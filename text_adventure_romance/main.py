"""Main entry point for the text adventure romance game."""

import sys

from text_adventure_romance.engine import GameState, display_text, run_chapter
from text_adventure_romance.story import (
    chapter1, chapter2, chapter3, chapter4, chapter5,
    ending_perfect, ending_good, ending_bittersweet, ending_sad,
)


def play():
    """Run the full game."""
    display_text("""
╔═══════════════════════════════════════════╗
║                                           ║
║     「 雨 天 的 邂 逅 」                  ║
║                                           ║
║     —— 一个关于遇见与心动的故事 ——        ║
║                                           ║
╚═══════════════════════════════════════════╝

  这是一个文字冒险恋爱游戏。
  你的每一个选择都会影响故事的走向和最终结局。
  请用心选择，用心感受。

  准备好了吗？按 Enter 开始……
""")

    try:
        input()
    except EOFError:
        pass

    state = GameState()

    chapters = [chapter1, chapter2, chapter3, chapter4, chapter5]
    for i, chapter_func in enumerate(chapters, 1):
        state.chapter = i
        run_chapter(chapter_func, state)

    # Determine and show ending
    ending_type = state.get_ending_type()
    endings = {
        "perfect": ending_perfect,
        "good": ending_good,
        "bittersweet": ending_bittersweet,
        "sad": ending_sad,
    }
    endings[ending_type](state)

    display_text("""
  感谢你玩完了这个游戏！
  不同的选择会带来不同的结局，
  欢迎再来一次，探索其他故事线。
""")


def main():
    """Entry point with error handling."""
    try:
        play()
    except KeyboardInterrupt:
        print("\n\n  再见，下次再来玩哦！\n")
        sys.exit(0)


if __name__ == "__main__":
    main()
