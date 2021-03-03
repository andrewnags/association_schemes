
.PHONY: check
check :
	mypy

.PHONY: fmt
fmt :
	flake8

.PHONY: test
test :
	pytest
