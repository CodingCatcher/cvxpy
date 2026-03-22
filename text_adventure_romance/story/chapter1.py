"""第一章：命运的邂逅 - 在雨中的咖啡馆遇见神秘的她。"""

from text_adventure_romance.engine import display_text, get_choice, get_player_name


def chapter1(state):
    display_text("""
═══════════════════════════════════════
  第一章：命运的邂逅
═══════════════════════════════════════

窗外下着淅淅沥沥的小雨。
你推开那家街角咖啡馆的门，铃铛发出清脆的响声。

咖啡馆里弥漫着浓郁的咖啡香气，
暖黄色的灯光让一切都变得柔和。
""")

    state.player_name = get_player_name()

    display_text(f"""
你走向吧台，准备点一杯热美式。
就在这时，你注意到角落里坐着一个女孩——
她低着头专注地在笔记本上画着什么，
一缕头发从耳边滑落，她浑然不觉。

她面前的咖啡似乎已经凉了很久。
""")

    display_text("你决定……")

    choice = get_choice([
        "主动走过去搭话",
        "默默坐在她旁边的位子",
        "帮她点一杯新的热咖啡送过去",
    ])

    state.add_choice(1, choice)

    if choice == 1:
        display_text("""
你鼓起勇气走到她面前。
"你好，我看你画得很专注……打扰一下可以吗？"

她抬起头，一双清澈的眼睛看着你，
先是有些惊讶，然后露出一个浅浅的微笑。

"啊……你好。没关系的，我只是在随便画画。"
她合上了笔记本，像是有些害羞。
""")
        state.modify_affection(10)
        state.set_flag("first_meeting_direct")

    elif choice == 2:
        display_text("""
你端着咖啡，悄悄坐在了她旁边的位子上。
你假装看着窗外的雨，余光却忍不住瞟向她的画。

过了一会儿，她似乎感觉到了什么，抬起头来。
你们的目光在空中相遇——

"……你一直在看我画画吗？"她歪了歪头。
你有些尴尬，但她却笑了起来。

"没关系啦，只是没人看过我画画，有点不习惯。"
""")
        state.modify_affection(5)
        state.modify_trust(10)
        state.set_flag("first_meeting_quiet")

    else:
        display_text("""
你走到吧台，多点了一杯热拿铁。
然后你端着两杯咖啡走向她的桌子。

"你的咖啡好像凉了，这杯请你。"

她惊讶地抬起头，眼中闪过一丝感动。
"谢谢你……你怎么知道我喜欢拿铁？"

"直觉吧。"你笑着说。

她低下头，嘴角微微上扬。
"谢谢，真的很暖。"
""")
        state.modify_affection(15)
        state.modify_trust(5)
        state.set_flag("first_meeting_coffee")

    display_text("""
你们聊了起来。
她叫林悦，是一个自由插画师。
她说她经常来这家咖啡馆寻找灵感。

"雨天的时候，世界好像会安静下来。"
她望着窗外说，"这种时候最适合画画了。"

时间过得很快，雨渐渐停了。
""")

    choice2 = get_choice([
        "要她的联系方式",
        '说"希望下次还能在这里遇见你"',
    ])

    state.add_choice(1, choice2)

    if choice2 == 1:
        display_text("""
"能加个微信吗？我很想看看你的其他作品。"

她犹豫了一秒，然后拿出手机。
"好啊……不过我画得不怎么样啦。"

你们交换了联系方式。
走出咖啡馆的时候，雨后的空气格外清新。
你看了看手机里新添的联系人——"林悦"。
心跳似乎比平时快了一些。
""")
        state.modify_affection(5)
        state.set_flag("got_contact")
    else:
        display_text("""
"希望下次还能在这里遇见你。"

她愣了一下，然后笑了。
"好啊，我经常在这里的。"

你走出咖啡馆，回头看了一眼。
她又低下头开始画画了，
但你注意到她的嘴角一直带着笑。

也许缘分这种东西，真的存在吧。
""")
        state.modify_trust(10)
        state.set_flag("left_to_fate")

    display_text("—— 第一章 完 ——")
