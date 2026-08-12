import numpy as np
import pandas as pd
import lightgbm as lgbm
import pickle
from pathlib import Path
import os

class Classifier_PAM50:
    """
    Классификатор PAM50.
    Модели: Basal, Her2, LumB, Normal. Краевой случай — LumA.
    На вход принимает один образец.
    """
    _pwd = Path(__file__).resolve().parent
    _models_path = os.path.join(_pwd, "subtype_classifiers_v4.pickle")

    def __init__(self):
        with open(self._models_path, "rb") as f:
            self.models = pickle.load(f)

    def _transform(self, x: pd.Series, features: list) -> pd.Series:
        """
        Возвращает образец, отмасштабированный по медиане и MAD.
        """
        x = x.loc[features].astype(float)
        m = x.median()
        mad = (x - m).abs().median()
        x_scaled = ((x - m) / mad).clip(-3, 3)
        return x_scaled

    def _predict_binary(self, subtype: str, sample: pd.Series):
        """
        Возвращает вывод бинарного классификатора
        """
        model_info = self.models[subtype]
        model = model_info["model"]
        feats = model_info["features"]

        x_tr = self._transform(sample, feats)
        X_tr = x_tr.to_frame().T

        proba = model.predict_proba(X_tr)[:, 1][0]
        y = model.predict(X_tr)[0]
        return y, proba

    def predict(self, sample: pd.Series) -> str:
        """
        Возвращает подтип PAM50
        """
        for subtype in ["Basal", "Her2", "LumB", "Normal"]:
            y, _ = self._predict_binary(subtype, sample)
            if y == 1:
                return subtype
        return "LumA"

    def predict_proba(self, sample: pd.Series) -> pd.Series:
        """
        Возвращает вектор вероятностей всех подтипов.
        """
        probs = {}
        for subtype in ["Basal", "Her2", "LumB", "Normal"]:
            _, p = self._predict_binary(subtype, sample)
            probs[subtype] = p

        probs["LumA"] = 1.0 - probs["Normal"]
        return pd.Series(probs, index=["Basal", "Her2", "LumB", "Normal", "LumA"])
