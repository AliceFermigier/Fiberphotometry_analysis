import matplotlib

def with_qt5agg():
    matplotlib.use("Qt5Agg", force=True)
    import matplotlib.pyplot as plt
    return plt

def with_agg():
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    return plt