## How to test the software?

We employ python [Unittest](https://docs.python.org/3/library/unittest.html) class to test our code. 
The file *test_suite.py* contains some Unittest Test Cases that should be run in order to be sure that the compilation went succesfully.

```bash
python3 test_suite.py
```

Or, for verbose output:
```bash
python3 test_suite.py -v
```

If everything went smoothly the output should look like the following:

```bash
...
Ran 20 tests in 0.140s

OK
```

*test_suite.py* makes use of the following packages:
* unittest
* pathlib
* os
* subprocess
