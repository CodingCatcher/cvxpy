"""Rich-based terminal UI for the game."""
from __future__ import annotations

import os
from typing import TYPE_CHECKING

from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.table import Table
from rich.columns import Columns
from rich import box

if TYPE_CHECKING:
    from .state import GameState
    from .characters import Character

console = Console()

CHAPTER_TITLES = {
    1: "第一章 · 最差劲的勇者",
    2: "第二章 · 低语之森",
    3: "第三章 · 圣城祭典",
    4: "第四章 · 魔族之女",
    5: "第五章 · 最终的选择",
}

HEART_FULL = "[red]\u2665[/red]"
HEART_EMPTY = "[dim]\u2661[/dim]"


def clear_screen() -> None:
    os.system("cls" if os.name == "nt" else "clear")


def show_title_screen() -> None:
    clear_screen()
    title = Text()
    title.append("\n\n")
    title.append("  \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\n", style="bright_yellow")
    title.append("  \u2588                                          \u2588\n", style="bright_yellow")
    title.append("  \u2588    \u4f20\u8bf4\u4e2d\u7684\u8272\u9b3c\u52c7\u8005\u7684\u6551\u8d4e          \u2588\n", style="bright_yellow")
    title.append("  \u2588    The Legendary Pervert's Redemption     \u2588\n", style="bright_yellow")
    title.append("  \u2588                                          \u2588\n", style="bright_yellow")
    title.append("  \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2588\n", style="bright_yellow")
    console.print(title, justify="center")
    console.print()
    console.print("[dim italic]    \u201c\u6211\u53ea\u662f\u60f3\u770b\u6f02\u4eae\u59d0\u59d0\u2026\u2026\u4e16\u754c\u548c\u5e73\u4ec0\u4e48\u7684\uff0c\u987a\u4fbf\u800c\u5df2\u3002\u201d\u2014\u2014\u96f7\u514b\u65af[/dim italic]", justify="center")
    console.print()


def show_menu(options: list[str]) -> int:
    console.print()
    for i, opt in enumerate(options, 1):
        console.print(f"  [bright_cyan][{i}][/bright_cyan] {opt}")
    console.print()
    while True:
        try:
            choice = console.input("[bright_cyan]\u8bf7\u9009\u62e9 > [/bright_cyan]")
            idx = int(choice.strip()) - 1
            if 0 <= idx < len(options):
                return idx
        except (ValueError, EOFError):
            pass
        console.print("[red]\u65e0\u6548\u7684\u9009\u62e9\uff0c\u8bf7\u91cd\u65b0\u8f93\u5165\u3002[/red]")


def show_chapter_header(chapter: int, location: str) -> None:
    title = CHAPTER_TITLES.get(chapter, f"\u7b2c{chapter}\u7ae0")
    header = f"{title}  \u2502  {location}" if location else title
    console.print(Panel(header, style="bright_yellow", box=box.DOUBLE))


def show_dialogue(speaker_name: str, speaker_color: str, text: str) -> None:
    name_text = Text(f"\u3010{speaker_name}\u3011", style=f"bold {speaker_color}")
    console.print()
    console.print(name_text)
    for line in text.split("\n"):
        console.print(f"  {line}")


def show_narration(text: str) -> None:
    console.print()
    console.print(Panel(text, style="italic dim white", box=box.SIMPLE, padding=(0, 2)))


def show_choices(choices: list[tuple[str, bool]]) -> int:
    """Show choices. Each tuple is (text, is_available). Returns index of chosen."""
    console.print()
    available_indices = []
    for i, (text, available) in enumerate(choices):
        if available:
            console.print(f"  [bright_cyan][{len(available_indices) + 1}][/bright_cyan] {text}")
            available_indices.append(i)
        else:
            console.print(f"  [dim strikethrough]  {text} (\u6761\u4ef6\u672a\u6ee1\u8db3)[/dim strikethrough]")
    console.print()
    while True:
        try:
            raw = console.input("[bright_cyan]\u4f60\u7684\u9009\u62e9 > [/bright_cyan]")
            idx = int(raw.strip()) - 1
            if 0 <= idx < len(available_indices):
                return available_indices[idx]
        except (ValueError, EOFError):
            pass
        console.print("[red]\u65e0\u6548\u7684\u9009\u62e9\uff0c\u8bf7\u91cd\u65b0\u8f93\u5165\u3002[/red]")


def show_affection_change(char_name: str, char_color: str, amount: int, new_total: int) -> None:
    if amount > 0:
        console.print(f"  [dim]\u2764 [{char_color}]{char_name}[/{char_color}] \u597d\u611f\u5ea6 +{amount} (\u2192 {new_total})[/dim]")
    elif amount < 0:
        console.print(f"  [dim]\u2764 [{char_color}]{char_name}[/{char_color}] \u597d\u611f\u5ea6 {amount} (\u2192 {new_total})[/dim]")


def show_affection_panel(state: GameState) -> None:
    from .characters import CHARACTERS
    table = Table(title="\u597d\u611f\u5ea6", box=box.ROUNDED, show_lines=False)
    table.add_column("\u89d2\u8272", style="bold")
    table.add_column("\u597d\u611f\u5ea6", justify="center")

    for char_id, value in state.affection.items():
        char = CHARACTERS.get(char_id)
        if not char:
            continue
        hearts = value // 10
        bar = HEART_FULL * hearts + HEART_EMPTY * (10 - hearts)
        table.add_row(f"[{char.color}]{char.name}[/{char.color}]", f"{bar} {value}")

    console.print()
    console.print(table)


def wait_for_continue() -> None:
    console.input("\n[dim]\u6309 Enter \u7ee7\u7eed...[/dim]")


def show_save_success(slot: int) -> None:
    console.print(f"\n[green]\u2713 \u5df2\u4fdd\u5b58\u5230\u5b58\u6863\u4f4d {slot}[/green]")


def show_game_over(ending_name: str, description: str) -> None:
    clear_screen()
    console.print()
    console.print(Panel(
        f"[bold]{ending_name}[/bold]\n\n{description}",
        title="\u2500\u2500 \u7ed3\u5c40 \u2500\u2500",
        style="bright_yellow",
        box=box.DOUBLE,
        padding=(1, 4),
    ))
    console.print()
    console.print("[dim italic]\u611f\u8c22\u60a8\u7684\u6e38\u73a9\uff01[/dim italic]", justify="center")
    console.print()
