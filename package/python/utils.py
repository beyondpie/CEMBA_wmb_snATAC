from inspect import getsource

def printfn(fn) -> str:
    """print fn source code"""
    content = getsource(fn)
    print(content)
    return(content)

def concat_allen_name(allen_l2):
    """allen_l2 is usually a pandas Series.
    """
    import re
    r = re.sub("_$", "", allen_l2.replace("/", "-").replace(" ", "_"))
    return(r)
