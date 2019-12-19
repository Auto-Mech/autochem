""" autoread parameters
"""
import autoparse.pattern as app


class Pattern():
    """ re patterns """
    ATOM_SYMBOL = (
        app.LETTER +
        app.maybe(app.LETTER)
    )
    NUMERIC_VALUE = app.NUMBER
