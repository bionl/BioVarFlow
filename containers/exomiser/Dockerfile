# Exomiser CLI on a base that has a shell and gcloud.
#
# The upstream exomiser/exomiser-cli image is distroless (jib): no /bin/sh, no
# package manager, entrypoint is the Java class directly. Nextflow wraps every
# script block in bash, so that image cannot back a Nextflow process at all --
# not the heredocs, not the data sync, not even echo. It also has no gcloud.
#
# /app/jib-classpath-file uses absolute /app paths, so the application copies
# across intact. The JVM is copied too (rather than apt-installed) to keep the
# exact version Exomiser was built against; /usr/lib/jvm is copied wholesale
# because the directory inside is arch-suffixed.
FROM google/cloud-sdk:slim

COPY --from=exomiser/exomiser-cli:15.1.0 /usr/lib/jvm /usr/lib/jvm
COPY --from=exomiser/exomiser-cli:15.1.0 /app /app

RUN ln -sf "$(ls -d /usr/lib/jvm/*/bin/java | head -1)" /usr/local/bin/java \
 && printf '#!/bin/bash\nexec java -cp @/app/jib-classpath-file org.monarchinitiative.exomiser.cli.Main "$@"\n' \
      > /usr/local/bin/exomiser \
 && chmod +x /usr/local/bin/exomiser

ENV EXOMISER_DATA_DIRECTORY=/exomiser-data

# Verified in this image:
#   bash 5.2.37 | gcloud | java 25.0.2 | exomiser 15.1.0
#   `exomiser batch --dry-run batch.txt` parses, resolves EXOMISER_DATA_DIRECTORY
#   and picks up EXOMISER_HG38_DATA_VERSION / EXOMISER_PHENOTYPE_DATA_VERSION.
#
# Build for linux/amd64 -- GCP Batch is x86, and the JVM directory copied above
# is arch-specific:
#   docker buildx build --platform linux/amd64 \
#     -t <registry>/exomiser-cli:15.1.0 --push containers/exomiser
