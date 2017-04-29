"""
    Create a console bar indicating the amount of time that's passed
"""

import sys


def setup(toolbar_width=40):
    """ Initialize console with toolbar_width spaces between [ ]

    Args:
        toolbar_width: Number of spaces between the brackets

    Returns:
        toolbar_width [int]"""

    # toolbar_width = 40
    global n_bar, n_done, upchar

    n_bar = toolbar_width
    n_done = 0
    upchar = '-'

    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1))

    return toolbar_width


def update(char='-'):
    """Prints char to indicate an update on the toolbar

    Args:
        char: The character to print"""

    global n_bar, n_done, upchar

    upchar = char

    if n_done < n_bar:
        sys.stdout.write(char)
        sys.stdout.flush()
    n_done += 1


def done(mapped=False):
    """Carriage return and flush the console"""

    global n_bar, n_done, upchar

    if mapped:
        sys.stdout.write("\r")
        sys.stdout.write("[%s]" % (upchar * n_bar))
        sys.stdout.write("\n")
        sys.stdout.flush()
    else:
        while n_done < n_bar:
            sys.stdout.write(upchar)
            sys.stdout.flush()
            n_done += 1

        sys.stdout.write("]")
        sys.stdout.write("\n")
        sys.stdout.flush()
