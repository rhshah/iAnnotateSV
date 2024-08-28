# The code defines custom exception classes for handling errors related to breakpoints in genetic
# annotation, such as breakpoints not in autosome or allosome and intergenic breakpoints.
class Error(Exception):
    """Base class for other exceptions"""

    pass


class ChrError(Error):
    """Raise when the breakpoint is not in autosome or allosome, where annotation is not possible"""

    def __init__(self, bkp):
        super().__init__(
            f"Breakpoint {bkp} is in neither autosome nor allosome, and cannot be annotated."
        )


class IntergenicError(Error):
    """Raise when a breakpoint in intergenic region cannot be resolved"""

    def __init__(self, bkp):
        super().__init__(f"Intergenic breakpoint {bkp} cannot be resolved.")
