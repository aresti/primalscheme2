class NoSuitablePrimers(Exception):
    """Unable to find any passing primers within the region."""

    pass


class NoReverseWindow(Exception):
    """
    Given the constraints of the next_pair (in scheme repair),
    the reverse window is insufficiently large to accomodate any
    primers.
    """

    pass
