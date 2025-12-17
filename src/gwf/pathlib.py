import pathlib


class Path(type(pathlib.Path())):
    _use_scratch: bool = False

    def __str__(self) -> str:
        return str(self)

    def __fspath__(self) -> str:
        return str(self)

    def __repr__(self):
        return str(self)

    @property
    def use_scratch(self) -> bool:
        return self._use_scratch


class TempPath(Path):
    pass
