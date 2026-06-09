process LOAD_VARIANTS {

    tag "${meta.sample} | ${meta.assay}"
    label 'low_cpu'

    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    tuple val(meta), path(variants_tsv)

    output:
    val(meta), emit: loaded_sample

    script:
    def sampleId   = meta.sample
    def gcsPath    = "${params.outdir}/variants/${variants_tsv.name}"
    def serviceUrl = params.vaic_service_url
    def apiKey     = params.vaic_api_key
    log.info "[LOAD_VARIANTS] ${sampleId} -> POST ${serviceUrl}/variants-db/load-variants  fileUrl=${gcsPath}"
    """
    python3 << 'PYEOF'
import requests, sys

resp = requests.post(
    "${serviceUrl}/variants-db/load-variants",
    json={"fileUrl": "${gcsPath}", "sample_id": "${sampleId}"},
    headers={"x-api-key": "${apiKey}"},
    timeout=600,
)
print(f"Status: {resp.status_code}")
print(resp.text)
resp.raise_for_status()
print("Variants loaded for ${sampleId}")
PYEOF
    """

    stub:
    """
    echo "stub: would load variants for ${meta.sample}"
    """
}
