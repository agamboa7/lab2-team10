# --- 1. Imports ---
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import matthews_corrcoef, classification_report, confusion_matrix, roc_curve, auc, precision_recall_curve, average_precision_score
import matplotlib.pyplot as plt
import seaborn as sns
import random
from pathlib import Path

# --- 2. Configuration ---
class Config:
    """Configuration class for all hyperparameters and settings."""
    # Data settings
    FILE_PATH = "merged_dataset_with_seqs.tsv"
    VAL_SET_SIZE = 0.15  # Use 15% of the training data for validation
    MAX_SEQ_LEN = 70     # Truncate/pad sequences to this length

    # Model hyperparameters
    INPUT_SIZE = 20      # One-hot encoding dimension for 20 amino acids
    CONV_OUT_CHANNELS = 64
    CONV_KERNEL_SIZE = 17
    LSTM_HIDDEN_SIZE = 256
    LSTM_NUM_LAYERS = 2
    DROPOUT_PROB = 0.5

    # Training settings
    BATCH_SIZE = 32
    LEARNING_RATE = 0.0001
    EPOCHS = 100
    PATIENCE = 15        # For early stopping
    MODEL_SAVE_PATH = "best_signal_peptide_model.pth"

    # Other settings
    RANDOM_SEED = 42
    DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# --- 3. Reproducibility ---
def set_seed(seed: int):
    """Sets the seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

# --- 4. Data Loading and Processing ---
def one_hot_encode(seq: str, max_len: int, aa_map: dict) -> np.ndarray:
    """One-hot encodes a single amino acid sequence with padding/truncation."""
    # Truncate sequence
    seq = seq[:max_len]

    # Create one-hot matrix
    one_hot = np.zeros((max_len, len(aa_map)), dtype=np.float32)
    for i, aa in enumerate(seq):
        if aa in aa_map:
            one_hot[i, aa_map[aa]] = 1.0

    return one_hot

class SignalPeptideDataset(Dataset):
    """Custom PyTorch Dataset for signal peptide sequences."""
    def __init__(self, sequences, labels):
        self.sequences = sequences
        self.labels = torch.tensor(labels, dtype=torch.float32)

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        return torch.tensor(self.sequences[idx], dtype=torch.float32), self.labels[idx]

def get_dataloaders(config: Config):
    """Loads data, creates splits, encodes, and returns DataLoaders."""
    print("--- Preparing Data ---")
    df = pd.read_csv(config.FILE_PATH, sep="\t")

    # Define the amino acid alphabet and mapping
    aa_alph = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aa_map = {aa: i for i, aa in enumerate(aa_alph)}

    # 1. Separate the untouchable test set
    df_test = df[df['Type'] == 'test']
    test_seqs = [one_hot_encode(s, config.MAX_SEQ_LEN, aa_map) for s in df_test['Sequence']]
    y_test = df_test['label'].tolist()

    # 2. Create a new validation set from the training data
    df_train_full = df[df['Type'] == 'train']
    train_full_seqs = df_train_full['Sequence'].tolist()
    train_full_labels = df_train_full['label'].tolist()

    X_train_seqs, X_val_seqs, y_train, y_val = train_test_split(
        train_full_seqs, train_full_labels,
        test_size=config.VAL_SET_SIZE,
        random_state=config.RANDOM_SEED,
        stratify=train_full_labels
    )

    # 3. One-hot encode the new train/val splits
    X_train = [one_hot_encode(s, config.MAX_SEQ_LEN, aa_map) for s in X_train_seqs]
    X_val = [one_hot_encode(s, config.MAX_SEQ_LEN, aa_map) for s in X_val_seqs]

    # 4. Create Datasets and DataLoaders
    train_dataset = SignalPeptideDataset(X_train, y_train)
    val_dataset = SignalPeptideDataset(X_val, y_val)
    test_dataset = SignalPeptideDataset(test_seqs, y_test)

    train_loader = DataLoader(train_dataset, batch_size=config.BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=config.BATCH_SIZE, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=config.BATCH_SIZE, shuffle=False)

    print(f"Train samples: {len(train_dataset)}, Val samples: {len(val_dataset)}, Test samples: {len(test_dataset)}")

    # 5. Calculate positive class weight for the loss function (from training set only)
    num_negatives = len(y_train) - sum(y_train)
    num_positives = sum(y_train)
    pos_weight = torch.tensor([num_negatives / num_positives], device=config.DEVICE)
    print(f"Calculated positive weight for loss function: {pos_weight.item():.2f}")

    return train_loader, val_loader, test_loader, pos_weight

# --- 5. Model Architecture ---
class SignalPeptideClassifier(nn.Module):
    """CNN + BiLSTM model for sequence classification."""
    def __init__(self, config: Config):
        super().__init__()
        self.conv_block = nn.Sequential(
            nn.Conv1d(
                in_channels=config.INPUT_SIZE,
                out_channels=config.CONV_OUT_CHANNELS,
                kernel_size=config.CONV_KERNEL_SIZE,
                padding='same'
            ),
            nn.BatchNorm1d(config.CONV_OUT_CHANNELS),
            nn.ReLU()
        )

        self.lstm = nn.LSTM(
            input_size=config.CONV_OUT_CHANNELS,
            hidden_size=config.LSTM_HIDDEN_SIZE,
            num_layers=config.LSTM_NUM_LAYERS,
            batch_first=True,
            bidirectional=True,
            dropout=config.DROPOUT_PROB if config.LSTM_NUM_LAYERS > 1 else 0
        )

        self.fc_block = nn.Sequential(
            nn.Linear(config.LSTM_HIDDEN_SIZE * 2, 512), # *2 for bidirectional
            nn.ReLU(),
            nn.Dropout(config.DROPOUT_PROB),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(config.DROPOUT_PROB),
            nn.Linear(256, 1) # Output logit
        )

    def forward(self, x):
        # x shape: (batch_size, seq_len, input_size)
        x = x.permute(0, 2, 1)  # (batch_size, input_size, seq_len) for Conv1d
        x = self.conv_block(x)
        x = x.permute(0, 2, 1)  # (batch_size, seq_len, conv_out) for LSTM

        out, _ = self.lstm(x)
        # We take the output of the last time step
        out = out[:, -1, :]

        logits = self.fc_block(out)
        return logits.squeeze(-1)

# --- 6. Trainer Class ---
class Trainer:
    """Handles the training, validation, and testing loops."""
    def __init__(self, model, train_loader, val_loader, test_loader, config, pos_weight):
        self.model = model.to(config.DEVICE)
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.test_loader = test_loader
        self.config = config
        self.device = config.DEVICE

        self.criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        self.optimizer = optim.Adam(self.model.parameters(), lr=config.LEARNING_RATE)

        self.best_val_mcc = -1.0
        self.epochs_no_improve = 0

    def _get_predictions(self, logits):
        """Converts logits to probabilities and binary predictions."""
        probs = torch.sigmoid(logits)
        preds = (probs > 0.5).float()
        return probs.cpu().numpy(), preds.cpu().numpy()

    def train(self):
        """Main training loop with validation and early stopping."""
        print("\n--- Starting Training ---")
        for epoch in range(self.config.EPOCHS):
            self.model.train()
            total_loss = 0
            for sequences, labels in self.train_loader:
                sequences, labels = sequences.to(self.device), labels.to(self.device)

                self.optimizer.zero_grad()
                logits = self.model(sequences)
                loss = self.criterion(logits, labels)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                self.optimizer.step()

                total_loss += loss.item()

            avg_train_loss = total_loss / len(self.train_loader)

            # Validation
            val_loss, val_mcc = self.evaluate(self.val_loader)

            print(f"Epoch [{epoch+1:03}/{self.config.EPOCHS:03}] | "
                  f"Train Loss: {avg_train_loss:.4f} | "
                  f"Val Loss: {val_loss:.4f} | "
                  f"Val MCC: {val_mcc:.4f}")

            # Early stopping and model checkpointing
            if val_mcc > self.best_val_mcc:
                self.best_val_mcc = val_mcc
                self.epochs_no_improve = 0
                torch.save(self.model.state_dict(), self.config.MODEL_SAVE_PATH)
                print(f"  -> Validation MCC improved. Saving model to '{self.config.MODEL_SAVE_PATH}'")
            else:
                self.epochs_no_improve += 1

            if self.epochs_no_improve >= self.config.PATIENCE:
                print(f"\nEarly stopping triggered after {self.config.PATIENCE} epochs without improvement.")
                break

    def evaluate(self, data_loader):
        """Evaluates the model on a given data loader."""
        self.model.eval()
        total_loss = 0
        all_labels = []
        all_preds = []
        with torch.no_grad():
            for sequences, labels in data_loader:
                sequences, labels = sequences.to(self.device), labels.to(self.device)

                logits = self.model(sequences)
                loss = self.criterion(logits, labels)
                total_loss += loss.item()

                _, preds = self._get_predictions(logits)

                all_labels.extend(labels.cpu().numpy())
                all_preds.extend(preds)

        avg_loss = total_loss / len(data_loader)
        mcc = matthews_corrcoef(all_labels, all_preds)
        return avg_loss, mcc

    def test(self):
        """Loads the best model and evaluates it on the test set."""
        print("\n--- Final Evaluation on Test Set ---")
        # Load the best model weights
        self.model.load_state_dict(torch.load(self.config.MODEL_SAVE_PATH))
        self.model.eval()

        all_labels, all_preds, all_probs = [], [], []
        with torch.no_grad():
            for sequences, labels in self.test_loader:
                sequences = sequences.to(self.device)
                logits = self.model(sequences)
                probs, preds = self._get_predictions(logits)

                all_labels.extend(labels.numpy())
                all_preds.extend(preds)
                all_probs.extend(probs)

        return np.array(all_labels), np.array(all_preds), np.array(all_probs)

# --- 7. Plotting and Reporting ---
def generate_report(y_true, y_pred, y_prob):
    """Generates a full classification report and visual plots."""
    print("\n" + "="*50)
    print("Classification Report")
    print("="*50)
    print(classification_report(y_true, y_pred, target_names=['No Signal Peptide', 'Signal Peptide']))

    # Create plots directory
    plots_dir = Path("plots")
    plots_dir.mkdir(exist_ok=True)

    # Confusion Matrix
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=['Predicted Negative', 'Predicted Positive'],
                yticklabels=['Actual Negative', 'Actual Positive'])
    plt.title('Confusion Matrix')
    plt.savefig(plots_dir / "confusion_matrix.png")
    plt.show()

    # ROC Curve
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    roc_auc = auc(fpr, tpr)
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    plt.savefig(plots_dir / "roc_curve.png")
    plt.show()

    # Precision-Recall Curve
    precision, recall, _ = precision_recall_curve(y_true, y_prob)
    pr_auc = average_precision_score(y_true, y_prob)
    plt.figure(figsize=(8, 6))
    plt.plot(recall, precision, color='blue', lw=2, label=f'PR curve (AP = {pr_auc:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend(loc="upper right")
    plt.grid(alpha=0.3)
    plt.savefig(plots_dir / "precision_recall_curve.png")
    plt.show()


if __name__ == "__main__":
    # Initialize configuration and set seed
    config = Config()
    set_seed(config.RANDOM_SEED)
    print(f"Using device: {config.DEVICE}")

    # Get data loaders and class weight
    train_loader, val_loader, test_loader, pos_weight = get_dataloaders(config)

    # Initialize model
    model = SignalPeptideClassifier(config)

    # Initialize and run trainer
    trainer = Trainer(model, train_loader, val_loader, test_loader, config, pos_weight)
    trainer.train()

    # Perform final testing and generate report
    y_true, y_pred, y_prob = trainer.test()
    generate_report(y_true, y_pred, y_prob)