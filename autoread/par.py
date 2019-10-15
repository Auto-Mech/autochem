""" autoread parameters
"""
import autoparse.pattern as app


class Pattern():
    """ re patterns """
    ATOM_SYMBOL = (
        app.LETTER +
        app.maybe(app.LETTER) # +
        # app.maybe(app.one_or_more(app.NUMBER))
    )
    NUMERIC_VALUE = app.NUMBER
