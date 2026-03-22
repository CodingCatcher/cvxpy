"""Text adventure romance game engine."""

import sys


class GameState:
    """Tracks player choices and relationship stats throughout the game."""

    def __init__(self):
        self.player_name = ""
        self.chapter = 1
        self.affection = 50  # 0-100 scale
        self.trust = 50
        self.choices_made = []
        self.flags = {}

    def add_choice(self, chapter, choice_id):
        self.choices_made.append((chapter, choice_id))

    def modify_affection(self, amount):
        self.affection = max(0, min(100, self.affection + amount))

    def modify_trust(self, amount):
        self.trust = max(0, min(100, self.trust + amount))

    def set_flag(self, flag_name, value=True):
        self.flags[flag_name] = value

    def get_flag(self, flag_name, default=False):
        return self.flags.get(flag_name, default)

    def get_ending_type(self):
        """Determine ending based on accumulated stats."""
        score = self.affection + self.trust
        if score >= 160:
            return "perfect"
        elif score >= 120:
            return "good"
        elif score >= 80:
            return "bittersweet"
        else:
            return "sad"


def display_text(text):
    """Display narrative text to the player."""
    print()
    for line in text.strip().split("\n"):
        print(f"  {line}")
    print()


def get_choice(options):
    """Present choices and get player input. Returns the choice index (1-based)."""
    for i, option in enumerate(options, 1):
        print(f"  [{i}] {option}")
    print()

    while True:
        try:
            raw = input("  你的选择 > ").strip()
            choice = int(raw)
            if 1 <= choice <= len(options):
                return choice
        except (ValueError, EOFError):
            pass
        print(f"  请输入 1-{len(options)} 之间的数字。")


def get_player_name():
    """Get the player's name."""
    while True:
        try:
            name = input("  请输入你的名字 > ").strip()
        except EOFError:
            name = ""
        if name:
            return name
        print("  名字不能为空哦。")


def run_chapter(chapter_func, state):
    """Run a single chapter function with the game state."""
    chapter_func(state)
