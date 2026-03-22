"""Tests for the game engine."""

import pytest
from text_adventure_romance.engine import GameState, display_text


class TestGameState:
    def test_initial_values(self):
        state = GameState()
        assert state.player_name == ""
        assert state.chapter == 1
        assert state.affection == 50
        assert state.trust == 50
        assert state.choices_made == []
        assert state.flags == {}

    def test_add_choice(self):
        state = GameState()
        state.add_choice(1, 2)
        state.add_choice(2, 1)
        assert state.choices_made == [(1, 2), (2, 1)]

    def test_modify_affection_normal(self):
        state = GameState()
        state.modify_affection(10)
        assert state.affection == 60

    def test_modify_affection_clamp_high(self):
        state = GameState()
        state.modify_affection(200)
        assert state.affection == 100

    def test_modify_affection_clamp_low(self):
        state = GameState()
        state.modify_affection(-200)
        assert state.affection == 0

    def test_modify_trust_normal(self):
        state = GameState()
        state.modify_trust(-10)
        assert state.trust == 40

    def test_modify_trust_clamp(self):
        state = GameState()
        state.modify_trust(100)
        assert state.trust == 100

    def test_flags(self):
        state = GameState()
        assert state.get_flag("test") is False
        state.set_flag("test", True)
        assert state.get_flag("test") is True
        state.set_flag("value", 42)
        assert state.get_flag("value") == 42

    def test_ending_perfect(self):
        state = GameState()
        state.affection = 90
        state.trust = 80
        assert state.get_ending_type() == "perfect"

    def test_ending_good(self):
        state = GameState()
        state.affection = 70
        state.trust = 60
        assert state.get_ending_type() == "good"

    def test_ending_bittersweet(self):
        state = GameState()
        state.affection = 50
        state.trust = 40
        assert state.get_ending_type() == "bittersweet"

    def test_ending_sad(self):
        state = GameState()
        state.affection = 30
        state.trust = 30
        assert state.get_ending_type() == "sad"

    def test_display_text(self, capsys):
        display_text("Hello\nWorld")
        captured = capsys.readouterr()
        assert "Hello" in captured.out
        assert "World" in captured.out


class TestEndingThresholds:
    """Test boundary conditions for ending determination."""

    def test_score_160_is_perfect(self):
        state = GameState()
        state.affection = 80
        state.trust = 80
        assert state.get_ending_type() == "perfect"

    def test_score_159_is_good(self):
        state = GameState()
        state.affection = 80
        state.trust = 79
        assert state.get_ending_type() == "good"

    def test_score_120_is_good(self):
        state = GameState()
        state.affection = 60
        state.trust = 60
        assert state.get_ending_type() == "good"

    def test_score_119_is_bittersweet(self):
        state = GameState()
        state.affection = 60
        state.trust = 59
        assert state.get_ending_type() == "bittersweet"

    def test_score_80_is_bittersweet(self):
        state = GameState()
        state.affection = 40
        state.trust = 40
        assert state.get_ending_type() == "bittersweet"

    def test_score_79_is_sad(self):
        state = GameState()
        state.affection = 40
        state.trust = 39
        assert state.get_ending_type() == "sad"
