class Error(Exception):
    '''Base class for other exceptions'''
    pass

class chrMT(Error):
    '''Raise when the breakpoint is in ChrMT, where annotation is not possible'''
    def __init__(self):
        Exception.__init__(
            self, "Breakpoint is in ChrMT, and cannot be annotated."
        )
