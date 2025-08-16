from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from imblearn.over_sampling import SMOTE
import seaborn as sns

# 创建保存图片的文件夹
def create_figure_dir(base_path):
    figure_dir = os.path.join(base_path, 'figure')
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)
    return figure_dir

# 创建保存 JSON 文件的文件夹
def create_json_dir(base_path):
    json_dir = os.path.join(base_path, 'json')
    if not os.path.exists(json_dir):
        os.makedirs(json_dir)
    return json_dir

# 绘制损失曲线
def plot_loss_curve(fold, epoch_steps, figure_dir, smoothing_window=50):
    train_losses_df = pd.read_csv(f'train_losses_fold_{fold}.csv')
    validate_losses_df = pd.read_csv(f'validate_losses_fold_{fold}.csv')

    steps_train = train_losses_df['Step']
    train_losses = train_losses_df['Loss']
    steps_validate = validate_losses_df['Step']
    validate_losses = validate_losses_df['Loss']

    smoothed_train_losses = train_losses.rolling(window=smoothing_window).mean()
    smoothed_validate_losses = validate_losses.rolling(window=smoothing_window).mean()

    plt.figure(figsize=(10, 5))
    plt.plot(steps_train, smoothed_train_losses, label='Smoothed Training Loss')
    plt.plot(steps_validate, smoothed_validate_losses, label='Smoothed Validation Loss')

    # 标注 epoch 起始位置
    for epoch, epoch_start in enumerate(epoch_steps):
        plt.axvline(x=epoch_start, color='grey', linestyle='--', alpha=0.7)
        plt.text(epoch_start, plt.ylim()[0], f'Epoch {epoch + 1}', rotation=90, verticalalignment='bottom', horizontalalignment='right', fontsize=8)

    plt.xlabel('Steps')
    plt.ylabel('Loss')
    plt.title(f'Loss Curve for Fold {fold}')
    plt.legend()
    plt.grid(True)

    # 确保 fold 文件夹存在
    fold_dir = os.path.join(figure_dir, f"fold_{fold}")
    os.makedirs(fold_dir, exist_ok=True)

    # 保存图片到对应的 fold 文件夹
    plt.savefig(os.path.join(fold_dir, f"loss_curve_fold_{fold}.svg"))
    plt.close()  # Use plt.close() instead of plt.show() for automatic saving without displaying

# 绘制 F1 分数曲线
def plot_f1_scores(f1_scores, fold, figure_dir):
    epochs = range(1, len(f1_scores["train"]) + 1)

    plt.figure(figsize=(10, 5))
    plt.plot(epochs, f1_scores["train"], label='Training F1 Score')
    plt.plot(epochs, f1_scores["validate"], label='Validation F1 Score')

    plt.xlabel('Epochs')
    plt.ylabel('F1 Score')
    plt.title(f'F1 Score Curve for Fold {fold}')
    plt.legend()
    plt.grid(True)

    # 确保 fold 文件夹存在
    fold_dir = os.path.join(figure_dir, f"fold_{fold}")
    os.makedirs(fold_dir, exist_ok=True)

    # 保存图片到对应的 fold 文件夹
    plt.savefig(os.path.join(fold_dir, f"f1_score_curve_fold_{fold}.svg"))
    plt.close()  # Use plt.close() instead of plt.show() for automatic saving without displaying

# 绘制混淆矩阵
def plot_confusion_matrix(y_true, y_pred, epoch, figure_dir, fold):
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(10, 7))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title(f'Confusion Matrix for Epoch {epoch}')

    # 确保 fold 文件夹存在
    fold_dir = os.path.join(figure_dir, f"fold_{fold}")
    os.makedirs(fold_dir, exist_ok=True)

    # 保存图片到对应的 fold 文件夹
    plt.savefig(os.path.join(fold_dir, f"confusion_matrix_epoch_{epoch}.svg"))
    plt.close()  # Use plt.close() instead of plt.show() for automatic saving without displaying

# 绘制 ROC 曲线
def plot_roc_curve(y_true, y_scores, epoch, figure_dir, fold):
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Receiver Operating Characteristic for Epoch {epoch}')
    plt.legend(loc='lower right')

    # 确保 fold 文件夹存在
    fold_dir = os.path.join(figure_dir, f"fold_{fold}")
    os.makedirs(fold_dir, exist_ok=True)

    # 保存图片到对应的 fold 文件夹
    plt.savefig(os.path.join(fold_dir, f"roc_curve_epoch_{epoch}.svg"))
    plt.close()  # Use plt.close() instead of plt.show() for automatic saving without displaying



data_folder = Path("data-space")
df = pd.read_csv(data_folder / "meta_genomic_data.csv")
classifiers = {
    "Random Forest": RandomForestClassifier(random_state=42, n_estimators=200, class_weight="balanced", max_depth=4),
    "SVM": SVC(random_state=42),
    "Gradient Boosting": GradientBoostingClassifier(random_state=42, max_depth=4, subsample=0.8),
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42)
}

X = df.iloc[:, 2:]
y = df['label']


ros = SMOTE(random_state=42)
X_resampled, y_resampled = ros.fit_resample(X, y)

classifiers = {
    "Random Forest": RandomForestClassifier(random_state=42, n_estimators=200, class_weight="balanced", max_depth=4),
    "SVM": SVC(random_state=42),
    "Gradient Boosting": GradientBoostingClassifier(random_state=42, max_depth=4, subsample=0.8),
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42)
}

results = {}

for name, model in classifiers.items():
    print(f"Evaluating {name}...")
    results[name] = evaluate_model(model, X_resampled, y_resampled)

for name, metrics in results.items():
    print(f"\n{name} Results:")
    print(f"Accuracy: {metrics['accuracy']:.4f}")
    print(f"Precision: {metrics['precision']:.4f}")
    print(f"Recall: {metrics['recall']:.4f}")
    print(f"F1 Score: {metrics['f1_score']:.4f}")

    plt.figure(figsize=(8, 6))
    sns.heatmap(metrics["confusion_matrix"], annot=True, fmt='d', cmap='Blues', xticklabels=np.unique(y),
                yticklabels=np.unique(y))
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title(f'{name} Confusion Matrix')
    plt.show()




