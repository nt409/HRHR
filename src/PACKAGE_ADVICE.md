

# How to setup python packages

Lots of stuff online, e.g. https://realpython.com/pypi-publish-python-package/

Some youtube stuff too. Try CookieCutter for a simple way?


## Updating existing package

Make changes, update version number (setup.py) and add to changelog

Build:
```bash
python setup.py sdist
```

Upload:
```bash
twine upload dist/*
```



