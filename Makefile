.PHONY: install test clean lint format docs

install:
	pip install -e .
	pip install -r requirements.txt

test:
	pytest Numerical-Optimization/tests/
	pytest Game-Theory/tests/
	pytest Convex-Optimization/tests/

clean:
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name "*.pyd" -delete
	find . -type f -name ".coverage" -delete
	find . -type d -name "*.egg-info" -exec rm -r {} +
	find . -type d -name "*.egg" -exec rm -r {} +
	find . -type d -name ".pytest_cache" -exec rm -r {} +
	find . -type d -name ".coverage" -exec rm -r {} +

lint:
	flake8 Numerical-Optimization/src/
	flake8 Game-Theory/src/
	black --check Numerical-Optimization/src/
	black --check Game-Theory/src/

format:
	black Numerical-Optimization/src/
	black Game-Theory/src/

docs:
	sphinx-build -b html docs/source docs/build/html 
