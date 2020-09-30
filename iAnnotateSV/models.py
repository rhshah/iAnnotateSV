class Error(Exception):
    '''Base class for other exceptions'''
    pass

class ChrError(Error):
    '''Raise when the breakpoint is not in autosome or allosome, where annotation is not possible'''
    def __init__(self, bkp):
        Exception.__init__(
                self, "Breakpoint " + str(bkp) + " is in neither autosome nor allosome, and cannot be annotated."
        )

class IntergenicError(Error):
    '''Raise when a breakpoint in intergenic region cannot be resolved'''
    def __init__(self, bkp):
        Exception.__init__(
            self, "Intergenic breakpoint " + str(bkp) + " cannot be resolved."
        )
