"""Game state management with save/load support."""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path

SAVE_DIR = Path.home() / ".adventure_save"
HEROINE_IDS = ["seraphina", "luna", "mira", "lilith", "elena"]


@dataclass
class GameState:
    current_scene: str = "prologue_01"
    chapter: int = 1
    affection: dict[str, int] = field(default_factory=lambda: {h: 0 for h in HEROINE_IDS})
    flags: set[str] = field(default_factory=set)
    history: list[str] = field(default_factory=list)

    def add_affection(self, char_id: str, amount: int) -> int | None:
        if char_id in self.affection:
            self.affection[char_id] = max(0, min(100, self.affection[char_id] + amount))
            return self.affection[char_id]
        return None

    def set_flag(self, flag: str) -> None:
        self.flags.add(flag)

    def has_flag(self, flag: str) -> bool:
        return flag in self.flags

    def get_top_heroine(self) -> str:
        return max(self.affection, key=self.affection.get)

    def save(self, slot: int = 1) -> Path:
        SAVE_DIR.mkdir(parents=True, exist_ok=True)
        path = SAVE_DIR / f"save_{slot}.json"
        data = {
            "current_scene": self.current_scene,
            "chapter": self.chapter,
            "affection": self.affection,
            "flags": list(self.flags),
            "history": self.history[-50:],
        }
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2))
        return path

    @classmethod
    def load(cls, slot: int = 1) -> GameState | None:
        path = SAVE_DIR / f"save_{slot}.json"
        if not path.exists():
            return None
        data = json.loads(path.read_text())
        return cls(
            current_scene=data["current_scene"],
            chapter=data["chapter"],
            affection=data["affection"],
            flags=set(data["flags"]),
            history=data["history"],
        )

    @classmethod
    def list_saves(cls) -> list[int]:
        if not SAVE_DIR.exists():
            return []
        saves = []
        for f in SAVE_DIR.glob("save_*.json"):
            try:
                saves.append(int(f.stem.split("_")[1]))
            except (ValueError, IndexError):
                pass
        return sorted(saves)
