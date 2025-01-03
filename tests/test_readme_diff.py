import re


def clean_content(content):
    """
    Removes logo-specific sections and whitespace for comparison.
    """
    # Remove all whitespace and normalize line endings
    content = content.replace("\r\n", "\n").strip()

    # Remove the logo sections using regex
    content = re.sub(r'<div align="center">.*?</div>', "", content, flags=re.DOTALL)
    return content


def test_readme_consistency():
    """
    Ensure GitHub and PyPI READMEs are consistent, excluding logos.
    """
    with open("README.md", "r", encoding="utf-8") as github_readme:
        github_content = github_readme.read()

    with open("README-pypi.md", "r", encoding="utf-8") as pypi_readme:
        pypi_content = pypi_readme.read()

    # Clean both
    cleaned_github = clean_content(github_content)
    cleaned_pypi = clean_content(pypi_content)

    # Assert they are identical
    assert (
        cleaned_github == cleaned_pypi
    ), "GitHub and PyPI README contents differ beyond logo sections!"
