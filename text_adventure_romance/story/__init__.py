"""Story chapters for the text adventure romance game."""

from .chapter1 import chapter1
from .chapter2 import chapter2
from .chapter3 import chapter3
from .chapter4 import chapter4
from .chapter5 import chapter5
from .endings import ending_perfect, ending_good, ending_bittersweet, ending_sad

__all__ = [
    "chapter1", "chapter2", "chapter3", "chapter4", "chapter5",
    "ending_perfect", "ending_good", "ending_bittersweet", "ending_sad",
]
