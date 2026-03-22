"""Core data model for scenes, choices, and the scene registry."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable


@dataclass
class Choice:
    """A player choice within a scene."""
    text: str
    next_scene: str
    affection: dict[str, int] = field(default_factory=dict)
    flags: list[str] = field(default_factory=list)
    condition: str | None = None


@dataclass
class Scene:
    """A single scene in the story."""
    id: str
    text: str
    speaker: str | None = None
    choices: list[Choice] | None = None
    next_scene: str | None = None
    on_enter: Callable | None = None
    chapter: int = 1
    location: str = ""


# Global scene registry
_scenes: dict[str, Scene] = {}


def register(*scenes: Scene) -> None:
    for s in scenes:
        _scenes[s.id] = s


def get_scene(scene_id: str) -> Scene:
    return _scenes[scene_id]


def all_scenes() -> dict[str, Scene]:
    return _scenes
