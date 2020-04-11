#!/bin/bash
echo "Building HTML output for IGoR asciidoc documentation..."
asciidoctor --doctype=article --backend=html --destination-dir=./docs/ --out-file=index.html -linkcss -a stylesheet=golo.css -a stylesdir=../ ./docs/asciidoc/IGoR_documentation.adoc
echo "Done."
