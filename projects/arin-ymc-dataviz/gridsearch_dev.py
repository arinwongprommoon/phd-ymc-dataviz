#!/usr/bin/env python3

# Example of constructing a pipeline and then perform a grid search.
# I can probably take better inspiration from Evangelos.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from classifier.binary import (
    PredictProbaHistogram,
    get_predictions,
    get_predictproba,
)
from classifier.gridsearch import GridSearchHandler
from classifier.metrics import ROCHandler, StratifiedKFoldHandler
from classifier.postpre import Catch22Transformer
from pathlib import Path
from postprocessor.core.processes.butter import butter
from postprocessor.core.processes.catch22 import catch22
from postprocessor.grouper import NameGrouper

# from pynndescent.sparse import sparse_correct_alternative_hellinger
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from slicer import slice_interval, slice_strain
from umapper import umap_plot

# 1. Construct pipeline, with classifier
my_pipeline = Pipeline(
    [
        ("featurise", Catch22Transformer()),
        ("scaler", StandardScaler()),
        # ("classifier", SVC(probability=True)),
        ("classifier", RandomForestClassifier()),
    ]
)

# 2. Define hyperparameter ranges to optimise over
window_range = [45]
C_range = [1]  # np.logspace(-3, 3, 3)
# gamma_range = np.logspace(-3, 3, 3)

param_grid = [
    {
        # "detrend__window": window_range,
        # "classifier__kernel": ["rbf"],
        # "classifier__C": C_range,
        #        'classifier__gamma': gamma_range
    },
]

# 3. Load features and labels
if True:
    # LOAD DATA/FEATURES
    data_directory = Path(
        # "/home/jupyter-arin/data/old-20220808/20016_2021_07_09_flavin_by4741_zwf1egf_01_00/",
        "/home/arin/git/data/20016_2021_07_09_flavin_by4741_zwf1egf_01_00/",
    )
    # Group by strain
    grouper = NameGrouper(data_directory)
    # Choose signals
    # TODO: Refactor & combine with strain_report_dev.py and dependencies to remove repetition
    signal_flavin = grouper.concat_signal("extraction/Flavin_bgsub/np_max/mean")
    signal_flavin_sliced = slice_interval(
        slice_strain(signal_flavin, ["by4741"]), 25, 168
    )
    signal_flavin_processed = butter.as_function(
        signal_flavin_sliced, critical_freq=1 / 350
    )
    signal_flavin_processed = signal_flavin_processed.droplevel(level=0)

    # LOAD LABELS
    targets_arin = pd.read_csv(
        "20016_by4741_labels.csv",
        index_col=[0, 1, 2],
    )
    targets_arin.columns = ["target"]

    features_arin = signal_flavin_processed.loc[targets_arin.index]  # .to_numpy()
    # targets_arin = targets_arin.to_numpy()

# THIS STEP SHOULD BE OPTIONAL
# 4. Perform grid search, then redefine hyperparameters for optimised pipeline
if False:
    gridsearcher = GridSearchHandler(param_grid)
    optimised_pipeline = gridsearcher.fit(
        pipeline=my_pipeline, features=features, targets=targets
    )

# 5. Split train-test and fit model
train_size = 0.75
features = features_arin
targets = targets_arin
features_train, features_test, targets_train, targets_test = train_test_split(
    features,
    targets,
    train_size=train_size,
)
optimised_pipeline = my_pipeline
optimised_pipeline.fit(features_train, targets_train)

# 5.1 Recursive feature elimination, with cross-validation
if False:
    rf = RandomForestClassifier()
    catch22_feature_matrix = catch22.as_function(features)
    scoring = "precision"

    min_features_to_select = 1
    rfecv = RFECV(
        estimator=rf,
        step=1,
        cv=StratifiedKFold(2),
        scoring=scoring,
        min_features_to_select=min_features_to_select,
    )
    rfecv.fit(features, targets)

    fig_rfe, ax_rfe = plt.subplots()
    ax_rfe.set_xlabel("Number of features selected")
    ax_rfe.set_ylabel("Cross validation score, " + scoring)
    ax_rfe.plot(
        range(min_features_to_select, len(rfecv.grid_scores_) + min_features_to_select),
        rfecv.grid_scores_,
    )

# 6. Get results
predictions_dict = get_predictions(
    optimised_pipeline, features_test, pd.unique(targets_arin.target)
)
targets_proba_df = get_predictproba(optimised_pipeline, features_test)
plotter = PredictProbaHistogram(
    category="oscillatory", predictproba_df=targets_proba_df
)
fig_proba, ax_proba = plt.subplots()
plotter.plot(ax_proba)

# 7. Evaluate
kfold = StratifiedKFoldHandler(optimised_pipeline, features, targets, 5)
kfold.pretty_print()
fig_kfold, ax_kfold = plt.subplots()
kfold.barplot(ax_kfold)

roc = ROCHandler(targets_proba_df, targets_test)
fig_roc, ax_roc = plt.subplots()
roc.plot(ax_roc)
print(roc.auc)

# 7.1 If Random Forest, get feature importances
feature_importances = pd.DataFrame(
    optimised_pipeline["classifier"].feature_importances_,
    # index=features.columns,
    columns=["importance"],
)
feature_importances_sorted = feature_importances.sort_values(
    ["importance"],
    ascending=False,
)
fig_imp, ax_imp = plt.subplots()
ax_imp.bar(
    x=np.arange(len(feature_importances_sorted)),
    height=feature_importances_sorted.to_numpy().T[0],
    tick_label=feature_importances_sorted.index,
)
# feature_importances_sorted

# TODO: Below, assumes that the 'strain' level exists.
# Not necessarily true.

# 8. Feature scatter plots
if False:
    catch22_feature_matrix = catch22.as_function(features)
    catch22_feature_matrix = catch22_feature_matrix.reset_index(level=["strain"])
    catch22_feature_matrix = catch22_feature_matrix.reset_index()
    # fig_pairplot, ax_pairplot = plt.subplots()
    sns.pairplot(catch22_feature_matrix, hue="strain", corner=True, plot_kws={"s": 1})

    corr = catch22_feature_matrix.drop(
        columns=[
            "cellID",
            "SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1",
            "SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1",
        ]
    ).corr()
    sns.clustermap(corr.abs())

# 9. UMAP
if False:
    umap_plot(
        data=catch22.as_function(features),
        n_neighbors=4,
        min_dist=1,
        n_components=2,
    )
