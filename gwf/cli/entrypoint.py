import sys

import colorama

from ..exceptions import GWFError
from .parsing import App


def main():
    colorama.init()

    app = App()
    try:
        app.run(sys.argv[1:])
    except GWFError as e:
        print("[Error] {}".format(str(e)))
        sys.exit(1)
