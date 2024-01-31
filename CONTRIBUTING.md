#Contributing to TACA code

When contribution to this package please have the following things in mind:

**NOTE**: _Please make sure that there are no exisiting [issues](<(https://github.com/SciLifeLab/TACA/issues)>) relating to whatever you want to report._

####To contribute:

1. Create an issue describing the bug / suggestion / improvement / ... [here](https://github.com/SciLifeLab/TACA/issues).
2. Fork this repository to your GitHub account
3. Make the necessary changes / additions to your forked TACA repository
4. Please _make sure_ that you've documented your code and changes using [sphinx](http://sphinx.readthedocs.org/en/latest/tutorial.html) syntax, as the documentation will be automatically generated using this engine, and published to [ReadTheDocs](http://project-management.readthedocs.org/)
5. Update the version number in `TACA/__init__.py`
6. Pull Request and wait for the responsible reviewer to review and merge the code

#### Version numbering

When contributing to this repo, also increment the version number in `TACA/__init__.py` as part of your pull request. If you suspect other pull requests will be merged before yours, simply ask the responsible reviewer to suggest a number after their review. Versioning follows the typical `MAJOR.MINOR.REVISION` system. Version numbers are set as follows:

- Revision - Small sized updates such as addressing a logic error or changing what version of a dependency is required
- Minor - Medium sized updates such as adding a new functionality or optimizing an existing feature
- Major - Large sized updates that either add a big bundle of functionality, or optimization that required a large section of the repository to be refactored.

That's it! And thanks a lot!
