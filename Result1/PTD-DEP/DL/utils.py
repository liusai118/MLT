import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
    precision_recall_fscore_support,
    roc_curve,
    auc,
)
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, DataStructs
from catboost import CatBoostClassifier
from imblearn.under_sampling import RandomUnderSampler
import optuna
from tqdm.notebook import tqdm


def get_best_threshold(tprs, fprs, thresholds, verbose=True):
    best_thresholds = []
    for tpr, fpr, threshs in zip(tprs, fprs, thresholds):
        j_scores = tpr - fpr
        best_threshold_index = np.argmax(j_scores)
        best_threshold = threshs[best_threshold_index]
        best_thresholds.append(best_threshold)
        if verbose:
            print(f"Best threshold: {best_threshold}")
    return np.mean(best_thresholds)


def get_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return fingerprint
    else:
        print(f"Invalid SMILES: {smiles}")
        return None


def filter_similar_molecules(
    df, smiles_col="SMILES", threshold=0.8, only_majority_class=False, target_class=None
):
    """
    Filter a dataframe to remove similar molecules based on a fingerprint similarity threshold.
    """
    df = df.copy()
    df["fingerprint"] = df[smiles_col].apply(get_fingerprint)
    df = df.dropna(subset=["fingerprint"])

    if only_majority_class:
        if target_class is None:
            raise ValueError(
                "target_class must be specified when only_majority_class=True."
            )
        majority_class = df[target_class].mode()[0]
        majority_df = df[df[target_class] == majority_class].copy()
        minority_df = df[df[target_class] != majority_class].copy()

        fingerprints = majority_df["fingerprint"].tolist()
        to_keep = [0]
        for i in tqdm(range(1, len(fingerprints)), desc="Filtering majority class"):
            similarities = np.array(
                DataStructs.BulkTanimotoSimilarity(
                    fingerprints[i], [fingerprints[j] for j in to_keep]
                )
            )
            if np.all(similarities < threshold):
                to_keep.append(i)

        filtered_majority_df = majority_df.iloc[to_keep]
        filtered_df = pd.concat([filtered_majority_df, minority_df])
    else:
        fingerprints = df["fingerprint"].tolist()
        to_keep = [0]
        for i in tqdm(range(1, len(fingerprints)), desc="Filtering molecules"):
            similarities = np.array(
                DataStructs.BulkTanimotoSimilarity(
                    fingerprints[i], [fingerprints[j] for j in to_keep]
                )
            )
            if np.all(similarities < threshold):
                to_keep.append(i)

        filtered_df = df.iloc[to_keep]

    filtered_df = filtered_df.drop(columns=["fingerprint"])
    return filtered_df


def check_imbalance(df, column, threshold=0.6):
    """Check if a binary column is imbalanced."""
    if df[column].value_counts(normalize=True).max() > threshold:
        return True
    return False


class DescriptorCalculator:
    def __init__(self, desc_list=None):
        self.desc_list = desc_list or Descriptors.descList

    def calculate_descriptors(self, mol):
        descriptors = {}
        for name, function in self.desc_list:
            descriptors[name] = function(mol)
        return descriptors


def count_atoms(mol, atom_symbol):
    return len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == atom_symbol])


count_nitrogen_atoms = lambda mol: count_atoms(mol, "N")
count_oxygen_atoms = lambda mol: count_atoms(mol, "O")


default_descriptors = [
    ("ExactMolWt", Descriptors.ExactMolWt),
    ("MolLogP", Descriptors.MolLogP),
    ("TPSA", Descriptors.TPSA),
    ("NumHDonors", Descriptors.NumHDonors),
    ("NumHAcceptors", Descriptors.NumHAcceptors),
    ("NumRotatableBonds", Descriptors.NumRotatableBonds),
    ("FractionCSP3", Descriptors.FractionCSP3),
    ("NumAromaticRings", Descriptors.NumAromaticRings),
    ("MaxPartialCharge", Descriptors.MaxPartialCharge),
    ("MinPartialCharge", Descriptors.MinPartialCharge),
    ("NumNitrogen", count_nitrogen_atoms),
    ("NumOxygen", count_oxygen_atoms),
]


def append_smiles_descriptors_to_df(
    df, smiles_col="SMILES", desc_list=None, keep_smiles_col=False
):
    """Append descriptors calculated from SMILES to a dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing a column with SMILES.
    smiles_col : str, optional
        Name of the column containing SMILES, by default "SMILES".
    desc_list : list, optional
        List of descriptors, by default None.
    keep_smiles_col : bool, optional
        Whether to keep the column containing SMILES, by default False.

    Returns
    -------
    pandas.DataFrame
        Dataframe with descriptors appended.
    """

    calculator = DescriptorCalculator(desc_list)

    problematic_smiles = []
    descriptors_list = []
    valid_indices = []

    for idx, smiles in df[smiles_col].items():
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise ValueError("Invalid molecule")
            descriptors = calculator.calculate_descriptors(mol)
            descriptors_list.append(descriptors)
            valid_indices.append(idx)
        except Exception as e:
            print(f"Problematic SMILES at index {idx}: {smiles}, Error: {e}")
            problematic_smiles.append((idx, smiles))

    descriptors_df = pd.DataFrame(descriptors_list)
    df = df.loc[valid_indices].reset_index(drop=True)
    df_with_descriptors = pd.concat([df, descriptors_df.reset_index(drop=True)], axis=1)

    if not keep_smiles_col:
        df_with_descriptors.drop(columns=[smiles_col], inplace=True)

    return df_with_descriptors, problematic_smiles


def append_morgan_fingerprints_to_df(
    df, smiles_col="SMILES", radius=2, nBits=2048, keep_smiles_col=False
):
    """
    Append Morgan fingerprints calculated from SMILES to a dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing a column with SMILES.
    smiles_col : str, optional
        Name of the column containing SMILES, by default "SMILES".
    radius : int, optional
        Radius for Morgan fingerprint, by default 2.
    nBits : int, optional
        Number of bits for Morgan fingerprint, by default 2048.
    keep_smiles_col : bool, optional
        Whether to keep the column containing SMILES, by default False.

    Returns
    -------
    pandas.DataFrame
        Dataframe with fingerprints appended.
    list
        List of tuples containing index and SMILES for problematic entries.
    """
    problematic_smiles = []
    fingerprints = []

    for idx, smiles in df[smiles_col].items():
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise ValueError("Invalid molecule")
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
            fingerprints.append(list(fingerprint))
        except Exception as e:
            print(f"Problematic SMILES at index {idx}: {smiles}, Error: {e}")
            problematic_smiles.append((idx, smiles))

    fingerprints_df = pd.DataFrame(
        fingerprints, columns=[f"Bit_{i}" for i in range(nBits)]
    )
    df_with_fingerprints = pd.concat([df, fingerprints_df], axis=1)

    if not keep_smiles_col:
        df_with_fingerprints.drop(columns=[smiles_col], inplace=True)

    return df_with_fingerprints, problematic_smiles


def run_experiments(
    create_model_func,
    X,
    y,
    n_splits=5,
    n_repeats=5,
    resampler=None,
    threshold=0.5,
    verbose=True,
):
    accs, precs, recs, f1s, results = [], [], [], [], []
    tprs = []  # True Positive Rates
    fprs = []  # False Positive Rates
    aucs = []  # Area Under the Curve
    thresholds = []
    aggregate_cm = np.zeros((2, 2), dtype=int)

    rskf = RepeatedStratifiedKFold(
        n_splits=n_splits, n_repeats=n_repeats, random_state=42
    )

    for i, (train_index, test_index) in enumerate(rskf.split(X, y)):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        clf = create_model_func()

        if resampler:
            X_train, y_train = resampler.fit_resample(X_train, y_train)

        clf.fit(X_train, y_train)

        # Make predictions on the test data
        y_proba = clf.predict_proba(X_test)[:, 1]

        # Apply the decision threshold
        y_pred = (y_proba >= threshold).astype(int)

        # ROC Curve statistics
        fpr, tpr, thr = roc_curve(y_test, y_proba)
        tprs.append(tpr)
        fprs.append(fpr)
        thresholds.append(thr)
        aucs.append(auc(fpr, tpr))

        # Confusion Matrix
        cm = confusion_matrix(y_test, y_pred)
        aggregate_cm += cm

        # Compute the metrics
        acc = accuracy_score(y_test, y_pred)
        prec = precision_score(y_test, y_pred)
        rec = recall_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)

        # Append the metrics for this fold to the lists
        accs.append(acc)
        precs.append(prec)
        recs.append(rec)
        f1s.append(f1)
        results.append((y_test, y_pred))

    if verbose:
        # Print the average metrics
        print(f"Average Accuracy: {np.mean(accs)}")
        print(f"Average Precision: {np.mean(precs)}")
        print(f"Average Recall: {np.mean(recs)}")
        print(f"Average F1 Score: {np.mean(f1s)}")
        print(aggregate_cm)

    return tprs, fprs, aucs, thresholds, aggregate_cm, results


def process_invalid_rows(df, remove_invalid=False):
    """
    Process a DataFrame to handle rows with infinite or NaN values.

    :param df: A pandas DataFrame to process.
    :param remove_invalid: A boolean flag to indicate whether to remove invalid rows.
    :return: A tuple containing the processed DataFrame and the DataFrame of invalid rows.
    """
    # Identify numeric columns
    numeric_columns = df.select_dtypes(include=np.number).columns

    # Initialize mask for invalid rows
    invalid_rows_mask = pd.Series(False, index=df.index)

    # Check for infinite values in numeric columns
    if not numeric_columns.empty:
        inf_mask = df[numeric_columns].apply(np.isinf).any(axis=1)
        invalid_rows_mask |= inf_mask

    # Check for NaN values in all columns
    nan_mask = df.isna().any(axis=1)
    invalid_rows_mask |= nan_mask

    # Identify invalid rows
    invalid_rows = df[invalid_rows_mask]

    # Removing invalid rows if requested
    if remove_invalid:
        df = df[~invalid_rows_mask]

    return df, invalid_rows


def finetune(
    X,
    y,
    min_iter=30,
    max_iter=200,
    min_depth=4,
    max_depth=10,
    n_trials=30,
    auto_class_weights=None,
):
    def objective(trial):
        # Define the hyperparameter search space
        n_estimators = trial.suggest_int("n_estimators", min_iter, max_iter)
        depth = trial.suggest_int("depth", min_depth, max_depth)

        # Modify the model creation function to use suggested parameters
        def create_model():
            return CatBoostClassifier(
                n_estimators=n_estimators,
                depth=depth,
                random_state=42,
                silent=True,
                auto_class_weights=auto_class_weights,
            )

        *_, results = run_experiments(
            create_model, X, y, resampler=None, threshold=0.5, verbose=False
        )

        # Calculate the average F1 score
        avg_f1 = np.mean([f1_score(y_true, y_pred) for y_true, y_pred in results])

        # Optuna tries to minimize the objective, so return the negative F1 score
        return -avg_f1

    pruner = optuna.pruners.MedianPruner(
        n_startup_trials=5, n_warmup_steps=10, interval_steps=1
    )

    # Create a study object with the pruner and specify the optimization direction
    study = optuna.create_study(direction="minimize", pruner=pruner)

    # Start the optimization
    study.optimize(objective, n_trials=n_trials)

    return study.best_params


def finetune_new(
    X,
    y,
    model_class,
    param_distributions,
    n_trials=30,
    resampler=None,
    threshold=0.5,
    verbose=False,
):
    def objective(trial):
        # Create a dictionary to hold the suggested hyperparameters
        params = {}
        for param_name, param_values in param_distributions.items():
            if param_values["type"] == "int":
                params[param_name] = trial.suggest_int(
                    param_name, param_values["low"], param_values["high"]
                )
            elif param_values["type"] == "float":
                params[param_name] = trial.suggest_float(
                    param_name, param_values["low"], param_values["high"]
                )
            elif param_values["type"] == "categorical":
                params[param_name] = trial.suggest_categorical(
                    param_name, param_values["choices"]
                )
        try:
            model = model_class(**params, random_state=42, silent=True)
        except Exception:
            model = model_class(**params, random_state=42)
        scores = cross_val_score(model, X, y, cv=5, scoring="f1_macro")
        avg_f1 = np.mean(scores)
        return -avg_f1

    pruner = optuna.pruners.MedianPruner(
        n_startup_trials=5, n_warmup_steps=10, interval_steps=1
    )

    # Create a study object with the pruner and specify the optimization direction
    study = optuna.create_study(direction="minimize", pruner=pruner)

    # Start the optimization
    study.optimize(objective, n_trials=n_trials)

    return study.best_params


def save_to_json(tprs, fprs, aucs, aggregate_cm, filename):
    # Convert numpy arrays to lists for JSON serialization
    tprs_list = [list(tpr) for tpr in tprs]
    fprs_list = [list(fpr) for fpr in fprs]
    aucs_list = list(aucs)
    aggregate_cm_list = aggregate_cm.tolist()

    # Pack data into a dictionary
    data = {
        "tprs": tprs_list,
        "fprs": fprs_list,
        "aucs": aucs_list,
        "aggregated_confusion_matrix": aggregate_cm_list,
    }

    # Write to JSON file
    with open(filename, "w") as file:
        json.dump(data, file)


def plot_roc(tprs, fprs, aucs, subfold_legend=False):
    for i in range(len(tprs)):
        plt.plot(
            fprs[i],
            tprs[i],
            lw=2,
            alpha=0.3,
            label=f"ROC fold {i+1} (AUC = {aucs[i]:.2f})" if subfold_legend else None,
        )

    # Calculating mean TPRs and FPRs for the mean ROC curve
    mean_fpr = np.linspace(0, 1, 100)
    mean_tpr = np.mean(
        np.array([np.interp(mean_fpr, fprs[i], tprs[i]) for i in range(len(tprs))]),
        axis=0,
    )
    mean_auc = auc(mean_fpr, mean_tpr)

    # Plotting the mean ROC curve
    plt.plot(
        mean_fpr,
        mean_tpr,
        color="blue",
        lw=2,
        alpha=0.8,
        label=f"Mean ROC (AUC = {mean_auc:.2f})",
    )

    # Final plot formatting
    plt.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", alpha=0.8)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curves")
    plt.legend(loc="lower right")
    plt.show()
