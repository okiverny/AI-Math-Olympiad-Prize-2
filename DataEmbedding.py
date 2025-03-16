from sentence_transformers import SentenceTransformer
from typing import List

class MathDataEmbedding:
    def __init__(self, model_name: str):
        self.model_name = model_name
        self.encoder_model = SentenceTransformer(model_name)

    def encode(self, texts: List[str], show_progress_bar: bool = False) -> List[List[float]]:
        return self.encoder_model.encode(texts, show_progress_bar=show_progress_bar)