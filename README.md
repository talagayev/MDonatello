MDonatello
==============================
[//]: # (Badges)

| **Latest release** | [![Last release tag][badge_release]][url_latest_release] ![GitHub commits since latest release (by date) for a branch][badge_commits_since]  [![Documentation Status][badge_docs]][url_docs]|
| :----------------- | :------- |
| **Workflows**      | ![Linux_CI_CD](https://github.com/talagayev/MDonatello/actions/workflows/Linux_CI_CD.yml/badge.svg) ![MacOS_CI_CD](https://github.com/talagayev/MDonatello/actions/workflows/MacOS_CI_CD.yml/badge.svg) ![Windows_CI_CD](https://github.com/talagayev/MDonatello/actions/workflows/Windows_CI_CD.yml/badge.svg)|
| **Status**         | [![codecov][badge_codecov]][url_codecov] [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)|
| **Community**      | [![License: GPL v2][badge_license]][url_license]  [![Powered by MDAnalysis][badge_mda]][url_mda]|

[badge_actions]: https://github.com/talagayev/mdonatello/actions/workflows/gh-ci.yaml/badge.svg
[badge_codecov]: https://codecov.io/gh/talagayev/mdonatello/branch/main/graph/badge.svg
[badge_commits_since]: https://img.shields.io/github/commits-since/talagayev/mdonatello/latest
[badge_docs]: https://readthedocs.org/projects/mdonatello/badge/?version=latest
[badge_license]: https://img.shields.io/badge/License-GPLv2-blue.svg
[badge_mda]: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
[badge_release]: https://img.shields.io/github/release-pre/talagayev/mdonatello.svg
[url_actions]: https://github.com/talagayev/mdonatello/actions?query=branch%3Amain+workflow%3Agh-ci
[url_codecov]: https://codecov.io/gh/talagayev/mdonatello/branch/main
[url_docs]: https://mdonatello.readthedocs.io/en/latest/?badge=latest
[url_latest_release]: https://github.com/talagayev/mdonatello/releases
[url_license]: https://www.gnu.org/licenses/gpl-2.0
[url_mda]: https://www.mdanalysis.org

2D small molecule visualization for MDAnalysis

This Repository is part of the Google Summer of Code 2024 for the following project:
[2D Visualization for small molecules](https://summerofcode.withgoogle.com/programs/2024/projects/sfy3kuqc)

The following Release is the submission for the GSoC 2024 Evaluation:
[GSoC 2024: 2D visualization for small molecules Release](https://github.com/talagayev/MDonatello/releases/tag/0.0.1)

Here is the final report for the GSoC2024 Project:
[GSoC 2024 - 2D visualization for small molecules](https://talagayev.github.io/posts/2024/08/1/)

MDonatello is bound by a [Code of Conduct](https://github.com/talagayev/mdonatello/blob/main/CODE_OF_CONDUCT.md).

### Installation

To build MDonatello from source,
First clone the repository:

```
git clone https://github.com/talagayev/MDonatello
```

Create a virtual environment and activate it:

```
conda create --name mdonatello
conda activate mdonatello
```

Then go into the MDonatello folder:

```
cd MDonatello
```

Finally this package from source:

```
pip install -e .
```

#### Running MDonatello

To use the **mdonatello** package you need to run a jupyter notebook, thus run the command:


```
jupyter notebook
```

Now that you started a jupyter notebook create a notebook file and enter the following command to use **mdonatello**:

```
import MDAnalysis as mda
import mdonatello
from mdonatello import MoleculeVisualizer

u = mda.Universe("input.pdb")
ag = u.select_atoms("resname UNK")
visualizer = MoleculeVisualizer(ag, show_atom_indices=False, width=-1, height=-1)
```

For a more detailed use of **MDonatello** and an example of the output follow the instructions in this section:


[Running MDonatello](https://mdonatello.readthedocs.io/en/latest/getting_started.html#usage)

### Copyright

The MDonatello source code is hosted at https://github.com/talagayev/mdonatello
and is available under the GNU General Public License, version 2 (see the file [LICENSE](https://github.com/talagayev/mdonatello/blob/main/LICENSE)).

Copyright (c) 2024, Valerij Talagayev


#### Acknowledgements
 
Project based on the 
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.
Please cite [MDAnalysis](https://github.com/MDAnalysis/mdanalysis#citation) when using MDonatello in published work.
