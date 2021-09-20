import re


def long_form(ref: str, que: str):
    cslong = []
    append = cslong.append
    _cs: str = ''
    _previous: str = ''
    for _ref, _que in zip(list(ref), list(que)):
        # Match
        if _ref == _que and _previous == "M":
            _cs = _ref
        elif _ref == _que and not _previous == "M":
            _cs = "=" + _ref
            _previous = "M"
        # Deletion
        elif _que == "D" and _previous == "D":
            _cs = _ref.lower()
        elif _que == "D" and not _previous == "D":
            _cs = "-" + _ref.lower()
            _previous = "D"
        # Insertion
        elif _ref == "I" and _previous == "I":
            _cs = _que.lower()
        elif _ref == "I" and not _previous == "I":
            _cs = "+" + _que.lower()
            _previous = "I"
        # Substitution
        elif _ref != _que:
            _cs = "*" + _ref.lower() + _que.lower()
            _previous = "S"
        append(_cs)
    return "cs:Z:" + ''.join(cslong)


def short_form(cslong):
    cs = []
    append = cs.append
    cs_split = re.split("(=[A-Z]+)", cslong)
    for _cs in cs_split:
        if _cs.startswith("="):
            _cs = ":" + str(len(re.findall("[A-Z]", _cs)))
        append(_cs)
    return ''.join(cs)
