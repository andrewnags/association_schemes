.PHONY: update
update :
	pip install -e .

.PHONY: check
check :
	mypy

.PHONY: lint
lint :
	flake8

.PHONY: test
test :
	pytest tests/*.py
