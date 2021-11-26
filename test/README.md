# A minimalistic but possibly helpful tests suite

To execute all tests declared in this directory simply run

	`./ALLTESTS.sh`

Principles:

* Tests are stored in `*.test` files.
* Every test consists in a single UNIX command, that must start with `TEST=<test-name>`
* A test is succesful if the exit code is zero.

In addition to `TEST`, a couple of environment variables are used to write the tests.
These are documented and initialised in `./ALLTESTS.sh` and used in the `*.test` files.


