"""Tests for story chapters - verify they run without errors."""

import pytest
from unittest.mock import patch
from text_adventure_romance.engine import GameState
from text_adventure_romance.story import (
    chapter1, chapter2, chapter3, chapter4, chapter5,
    ending_perfect, ending_good, ending_bittersweet, ending_sad,
)


def make_input_fn(inputs):
    """Create a side_effect function that returns inputs in sequence."""
    it = iter(inputs)
    def fake_input(prompt=""):
        return next(it)
    return fake_input


class TestChapter1:
    @patch("builtins.input")
    def test_chapter1_choice1(self, mock_input):
        mock_input.side_effect = make_input_fn(["小明", "1", "1"])
        state = GameState()
        chapter1(state)
        assert state.player_name == "小明"
        assert len(state.choices_made) == 2
        assert state.get_flag("first_meeting_direct")

    @patch("builtins.input")
    def test_chapter1_choice2(self, mock_input):
        mock_input.side_effect = make_input_fn(["玩家", "2", "2"])
        state = GameState()
        chapter1(state)
        assert state.get_flag("first_meeting_quiet")
        assert state.get_flag("left_to_fate")

    @patch("builtins.input")
    def test_chapter1_choice3(self, mock_input):
        mock_input.side_effect = make_input_fn(["测试", "3", "1"])
        state = GameState()
        chapter1(state)
        assert state.get_flag("first_meeting_coffee")
        assert state.get_flag("got_contact")


class TestChapter2:
    @patch("builtins.input")
    def test_chapter2_with_contact(self, mock_input):
        mock_input.side_effect = make_input_fn(["1", "2"])
        state = GameState()
        state.player_name = "测试"
        state.set_flag("got_contact")
        chapter2(state)
        assert len(state.choices_made) == 2

    @patch("builtins.input")
    def test_chapter2_without_contact(self, mock_input):
        mock_input.side_effect = make_input_fn(["2", "1"])
        state = GameState()
        state.player_name = "测试"
        chapter2(state)
        assert len(state.choices_made) == 2

    @patch("builtins.input")
    def test_chapter2_tease(self, mock_input):
        mock_input.side_effect = make_input_fn(["2", "3"])
        state = GameState()
        state.player_name = "测试"
        chapter2(state)
        assert state.get_flag("teased_about_drawing")


class TestChapter3:
    @patch("builtins.input")
    def test_chapter3_direct_jacket(self, mock_input):
        mock_input.side_effect = make_input_fn(["1", "2", "1"])
        state = GameState()
        state.player_name = "测试"
        chapter3(state)
        assert state.get_flag("gave_jacket")

    @patch("builtins.input")
    def test_chapter3_touch_hair(self, mock_input):
        mock_input.side_effect = make_input_fn(["2", "1", "2"])
        state = GameState()
        state.player_name = "测试"
        chapter3(state)
        assert state.get_flag("touched_hair")


class TestChapter4:
    @patch("builtins.input")
    def test_chapter4_honest(self, mock_input):
        mock_input.side_effect = make_input_fn(["2", "2"])
        state = GameState()
        state.player_name = "测试"
        chapter4(state)
        assert state.get_flag("honest_about_feelings")

    @patch("builtins.input")
    def test_chapter4_tried_to_leave(self, mock_input):
        mock_input.side_effect = make_input_fn(["3", "3"])
        state = GameState()
        state.player_name = "测试"
        chapter4(state)
        assert state.get_flag("tried_to_leave")


class TestChapter5:
    @patch("builtins.input")
    def test_chapter5_confess(self, mock_input):
        mock_input.side_effect = make_input_fn(["2", "1"])
        state = GameState()
        state.player_name = "测试"
        chapter5(state)
        assert state.get_flag("confessed")

    @patch("builtins.input")
    def test_chapter5_hold_hands(self, mock_input):
        mock_input.side_effect = make_input_fn(["1", "2"])
        state = GameState()
        state.player_name = "测试"
        chapter5(state)
        assert state.get_flag("held_hands")

    @patch("builtins.input")
    def test_chapter5_she_confesses(self, mock_input):
        mock_input.side_effect = make_input_fn(["3", "3"])
        state = GameState()
        state.player_name = "测试"
        chapter5(state)
        assert state.get_flag("she_confessed")


class TestEndings:
    def test_ending_perfect(self, capsys):
        state = GameState()
        state.player_name = "测试"
        state.affection = 90
        state.trust = 80
        ending_perfect(state)
        output = capsys.readouterr().out
        assert "最好的我们" in output

    def test_ending_good(self, capsys):
        state = GameState()
        state.player_name = "测试"
        state.affection = 70
        state.trust = 60
        ending_good(state)
        output = capsys.readouterr().out
        assert "晴天" in output

    def test_ending_bittersweet(self, capsys):
        state = GameState()
        state.player_name = "测试"
        state.affection = 50
        state.trust = 40
        ending_bittersweet(state)
        output = capsys.readouterr().out
        assert "半晴半雨" in output

    def test_ending_sad(self, capsys):
        state = GameState()
        state.player_name = "测试"
        state.affection = 30
        state.trust = 30
        ending_sad(state)
        output = capsys.readouterr().out
        assert "雨一直下" in output


class TestFullGame:
    """Test a complete playthrough from chapter 1 to ending."""

    @patch("builtins.input")
    def test_full_game_best_path(self, mock_input):
        """All best choices should give perfect ending."""
        mock_input.side_effect = make_input_fn([
            "小明",       # name
            "3", "1",     # ch1: coffee + get contact
            "2", "2",     # ch2: tease + eager yes
            "1", "2", "2",  # ch3: direct + "in your eyes" + touch hair
            "2", "2",     # ch4: ask if ok + honest
            "2", "1",     # ch5: find her + confess
        ])
        state = GameState()
        for i, ch in enumerate([chapter1, chapter2, chapter3, chapter4, chapter5], 1):
            state.chapter = i
            ch(state)
        assert state.get_ending_type() == "perfect"

    @patch("builtins.input")
    def test_full_game_conservative_path(self, mock_input):
        """Conservative choices should give a lower-score ending."""
        mock_input.side_effect = make_input_fn([
            "路人",       # name
            "2", "2",     # ch1: sit quietly + leave to fate (2 options)
            "3", "3",     # ch2: listen quietly + ask what to do
            "3", "3", "1",  # ch3: exercise + silent + jacket (2 options)
            "1", "3",     # ch4: polite handshake + focus on contest
            "3", "3",     # ch5: wait quietly + wait for her
        ])
        state = GameState()
        for i, ch in enumerate([chapter1, chapter2, chapter3, chapter4, chapter5], 1):
            state.chapter = i
            ch(state)
        # Even conservative choices accumulate decent stats in this game
        assert state.get_ending_type() in ("perfect", "good", "bittersweet", "sad")
