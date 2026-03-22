"""Core game engine: scene management, choice handling, game loop."""
from __future__ import annotations

from .state import GameState
from .story.common import get_scene, Scene, Choice
from .characters import CHARACTERS, NARRATOR, get_char
from . import ui


class GameEngine:
    def __init__(self, state: GameState | None = None):
        self.state = state or GameState()
        self.running = True

    def run(self) -> None:
        while self.running:
            scene = get_scene(self.state.current_scene)
            self._play_scene(scene)

    def _play_scene(self, scene: Scene) -> None:
        if scene.chapter and scene.chapter != self.state.chapter:
            self.state.chapter = scene.chapter
            ui.clear_screen()

        ui.show_chapter_header(self.state.chapter, scene.location)

        if scene.on_enter:
            scene.on_enter(self.state)

        self.state.history.append(scene.id)

        if scene.speaker:
            char = get_char(scene.speaker)
            ui.show_dialogue(char.name, char.color, scene.text)
        else:
            ui.show_narration(scene.text)

        if scene.choices:
            self._handle_choices(scene)
        elif scene.next_scene:
            self._handle_continue(scene)
        else:
            self.running = False
            return

    def _handle_choices(self, scene: Scene) -> None:
        choices = scene.choices
        choice_display = []
        for c in choices:
            available = True
            if c.condition and not self.state.has_flag(c.condition):
                available = False
            choice_display.append((c.text, available))

        # Add system options
        choice_display.append(("\u67e5\u770b\u597d\u611f\u5ea6", True))
        choice_display.append(("\u4fdd\u5b58\u6e38\u620f", True))

        while True:
            idx = ui.show_choices(choice_display)

            # System options
            if idx == len(choices):
                ui.show_affection_panel(self.state)
                ui.show_chapter_header(self.state.chapter, scene.location)
                if scene.speaker:
                    char = get_char(scene.speaker)
                    ui.show_dialogue(char.name, char.color, scene.text)
                else:
                    ui.show_narration(scene.text)
                continue
            elif idx == len(choices) + 1:
                self._save_menu()
                ui.show_chapter_header(self.state.chapter, scene.location)
                if scene.speaker:
                    char = get_char(scene.speaker)
                    ui.show_dialogue(char.name, char.color, scene.text)
                else:
                    ui.show_narration(scene.text)
                continue

            chosen = choices[idx]
            self._apply_choice(chosen)
            self.state.current_scene = chosen.next_scene
            break

    def _handle_continue(self, scene: Scene) -> None:
        ui.wait_for_continue()
        self.state.current_scene = scene.next_scene

    def _apply_choice(self, choice: Choice) -> None:
        for flag in choice.flags:
            self.state.set_flag(flag)

        for char_id, amount in choice.affection.items():
            new_val = self.state.add_affection(char_id, amount)
            if new_val is not None:
                char = get_char(char_id)
                ui.show_affection_change(char.name, char.color, amount, new_val)

    def _save_menu(self) -> None:
        ui.console.print("\n[bright_yellow]\u9009\u62e9\u5b58\u6863\u4f4d:[/bright_yellow]")
        slot = 1
        try:
            raw = ui.console.input("[bright_cyan]\u5b58\u6863\u4f4d (1-3) > [/bright_cyan]")
            slot = int(raw.strip())
            if slot < 1 or slot > 3:
                slot = 1
        except (ValueError, EOFError):
            pass
        path = self.state.save(slot)
        ui.show_save_success(slot)
