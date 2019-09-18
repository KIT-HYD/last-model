# LAST-model

## older MATLAB & Python version

The version used in the publication XXX, written in MATLAB, as well as its port to Python
are currently kept in the example folder for reference. This Model is under heavy development
and will substantially change in the future. The older code will, however, be kept, but 
maybe moved to another location. This will be stated here, in case you have been redirected here
from the aforementioned publication.

## LAST

### Install

You can install LAST using pip:

```bash
pip install last_model
```

A more recent (and maybe more unstable) version might be available on Github.

```bash
git clone https://github.com/KIT-HYD/last_model.git
cd last_model
pip install -r requirements.txt
python setup.py install
```

Using pip is the recommended way.

## Using LAST

We will prepare a documentation soon. The Last model can be instantiated 
using the `Last` class and executed using the `Last.run` method. 
Model setup and model parameters are specified in a JSON file. A sample 
test file is supplied with `last_model`.

```python
from last_model import Last

last = Last()
last.run()
```

## Notes:

* Last is slow at the moment. We are working on an improved performance
* Last is undocumented. This will be changed.
* Last is untested. We will provide unit tests with a future release
* Last is pre-alpha. Interfaces and internal workflow will change substantially.
* Last is not verbose, but it will learn speaking in a future release
