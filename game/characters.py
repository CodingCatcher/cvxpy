"""Character definitions for the game."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Character:
    id: str
    name: str
    title: str
    color: str
    description: str


CHARACTERS: dict[str, Character] = {}


def _reg(c: Character) -> Character:
    CHARACTERS[c.id] = c
    return c


REX = _reg(Character(
    id="rex",
    name="雷克斯 (Rex)",
    title="传说中的色鬼冒险者",
    color="bright_yellow",
    description="实力超群却好色成性的冒险者。被每一个队伍踢出来过。",
))

SERAPHINA = _reg(Character(
    id="seraphina",
    name="赛拉菲娜 (Seraphina)",
    title="圣骑士团副团长",
    color="gold1",
    description="冷傲高贵的女圣骑士，被指派监视雷克斯。外冷内热。",
))

LUNA = _reg(Character(
    id="luna",
    name="露娜 (Luna)",
    title="失忆的精灵法师",
    color="cyan",
    description="在低语之森中被发现的神秘精灵，失去了大部分记忆。",
))

MIRA = _reg(Character(
    id="mira",
    name="米拉 (Mira)",
    title="阳光牧师",
    color="bright_green",
    description="圣城神殿的见习牧师，天真开朗但洞察力惊人。",
))

LILITH = _reg(Character(
    id="lilith",
    name="莉莉丝 (Lilith)",
    title="半魔族叛逃者",
    color="magenta",
    description="从魔王军叛逃的半魔族少女，孤独而倔强。",
))

ELENA = _reg(Character(
    id="elena",
    name="艾莲娜 (Elena)",
    title="青梅竹马的酒馆老板娘",
    color="red",
    description="雷克斯的青梅竹马，经营着一家小酒馆，暗恋雷克斯多年。",
))

NARRATOR = Character(
    id="narrator",
    name="旁白",
    title="",
    color="white",
    description="",
)

GUILD_MASTER = Character(
    id="guild_master",
    name="公会长·巴尔德 (Bard)",
    title="冒险者公会长",
    color="bright_blue",
    description="冒险者公会的老会长，对雷克斯又爱又恨。",
)

DEMON_LORD = Character(
    id="demon_lord",
    name="魔王·阿扎泽尔 (Azazel)",
    title="黑暗之王",
    color="bright_red",
    description="统治魔族的强大存在，似乎与露娜有某种联系。",
)


def get_char(char_id: str) -> Character:
    return CHARACTERS.get(char_id, NARRATOR)
