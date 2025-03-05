import ast
import pathlib
import textwrap

import yaml

api_docs = pathlib.Path("docs/api")
api_docs.mkdir(exist_ok=True)

src = pathlib.Path("pypythia")

api_files = {
    "msa": None,
    "raxmlng": None,
    "prediction": None,
    "predictor": None,
    "custom_types": None,
    "custom_errors": None,
    "config": None,
}

for file_name in api_files:
    file = src / f"{file_name}.py"
    api_file = api_docs / (file_name + ".md")

    api_files[file_name] = api_file

    node = ast.parse(file.read_text())

    classes = []
    methods = []

    for item in node.body:
        if isinstance(item, ast.FunctionDef):
            fn_name = item.name
            if fn_name.startswith("_"):
                continue
            methods.append(item.name)
        elif isinstance(item, ast.ClassDef):
            classes.append(item.name)

    with api_file.open("w") as f:
        for cls in classes:
            f.write(
                textwrap.dedent(
                    f"""
            ::: pypythia.{file.stem}.{cls}\n
                options:
                    show_root_heading: true
                    merge_init_into_class: false
                    group_by_category: true
                    modernize_annotations: true
            """
                )
            )

        for mtd in methods:
            f.write(
                textwrap.dedent(
                    f"""
            ::: pypythia.{file.stem}.{mtd}\n
                options:
                    show_root_heading: true
                    modernize_annotations: true
            """
                )
            )

    if file_name == "config":
        with api_file.open("a") as f:
            f.write(
                textwrap.dedent(
                    f"""
            ::: pypythia.{file.stem}.DEFAULT_MODEL_FILE\n
                options:
                    show_root_heading: true
                    modernize_annotations: true

            ::: pypythia.{file.stem}.DEFAULT_RAXMLNG_EXE\n
                options:
                    show_root_heading: true
                    modernize_annotations: true
            """
                )
            )

mkdocs_cfg_file = pathlib.Path("mkdocs.yml")
mkdocs_cfg = yaml.safe_load(mkdocs_cfg_file.read_text())


api_nav = {"API Reference": [{name: f"api/{f.name}"} for name, f in api_files.items()]}

nav = []

for el in mkdocs_cfg["nav"]:
    if "API Reference" in el:
        continue
    nav.append(el)

nav.append(api_nav)
mkdocs_cfg["nav"] = nav

yaml.dump(mkdocs_cfg, mkdocs_cfg_file.open("w"))
