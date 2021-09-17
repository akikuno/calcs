import re

ref = "AIIIAAAAAATT"
que = "ATTTDDDAAAAA"


def call_cs_long(ref, que):
    cslong = []
    append = cslong.append
    _cs = ''
    _previous = ''
    for _ref, _que in zip(list(ref), list(que)):
        if _ref == _que:
            # Match
            if _previous == "M":
                _cs = _ref
            else:
                _cs = "=" + _ref
            _previous = "M"
        # Deletion
        elif _que == "D":
            if _previous == "D":
                _cs = _ref.lower()
            else:
                _cs = "-" + _ref.lower()
            _previous = "D"
        # Insertion
        elif _ref == "I":
            if _previous == "I":
                _cs = _que.lower()
            else:
                _cs = "+" + _que.lower()
            _previous = "I"
        # Substitution
        elif _ref != _que:
            _cs = "*" + _ref.lower() + _que.lower()
        append(_cs)
    return "cs:Z:" + ''.join(cslong)


def call_cs_short(cslong):
    cs = []
    append = cs.append
    cs_split = re.split("(=[A-Z]+)", cslong)
    for _cs in cs_split:
        if _cs.startswith("="):
            _cs = ":" + str(len(re.findall("[A-Z]", _cs)))
        append(_cs)
    return ''.join(cs)


cslong = call_cs_long(ref, que)

cs = call_cs_short(cslong)

cslong
cs
