def sort2lists(lst1, lst2):
    """
    sorts two unsorted lists according to first given list
    Parameters
    ----------
    lst1: list of sortables
    lst2: list

    Returns
    -------
    l1, l2: sorted lists

    """
    return [l1 for l1, _ in sorted(zip(lst1,lst2))], [l2 for _, l2 in sorted(zip(lst1,lst2))]

