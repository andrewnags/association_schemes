
.PHONY: check
check :
	mypy

.PHONY: lint
lint :
	flake8

.PHONY: test
test :
	pytest
