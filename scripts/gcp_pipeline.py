#!/usr/bin/env python3
"""Google Cloud Life Sciences pipeline for NRXN1 analysis."""

from google.cloud import storage
from google.cloud import lifesciences_v2beta as lifesciences
from typing import Dict, Any, Optional, List
from pathlib import Path
import time
import json


class GCPPipelineRunner:
    """Run NRXN1 analysis on Google Cloud Life Sciences."""
    
    def __init__(
        self,
        project_id: str,
        location: str = "us-central1",
        service_account: Optional[str] = None
    ):
        """Initialize GCP pipeline runner.
        
        Args:
            project_id: GCP project ID.
            location: GCP region.
            service_account: Optional service account email.
        """
        self.project_id = project_id
        self.location = location
        self.service_account = service_account
        
        self.lifesciences_client = lifesciences.LifeSciencesClient()
        self.storage_client = storage.Client(project=project_id)
        
        self.parent = f"projects/{project_id}/locations/{location}"
    
    def run_pipeline(
        self,
        input_gcs_path: str,
        output_gcs_path: str,
        machine_type: str = "n1-standard-8",
        disk_size_gb: int = 100,
        docker_image: str = "gcr.io/your-project/nrxn1-analysis:latest"
    ) -> str:
        """Run the analysis pipeline.
        
        Args:
            input_gcs_path: GCS path to input VCF.
            output_gcs_path: GCS path for output.
            machine_type: GCE machine type.
            disk_size_gb: Boot disk size.
            docker_image: Docker image to use.
        
        Returns:
            Operation name.
        """
        pipeline = lifesciences.Pipeline(
            actions=[
                lifesciences.Action(
                    image_uri="google/cloud-sdk:slim",
                    commands=[
                        "gsutil", "cp", input_gcs_path, "/data/input.vcf"
                    ],
                    mounts=[
                        lifesciences.Mount(disk="data", path="/data")
                    ]
                ),
                lifesciences.Action(
                    image_uri=docker_image,
                    commands=[
                        "python", "main.py",
                        "--input", "/data/input.vcf",
                        "--output", "/data/output",
                        "--mode", "full"
                    ],
                    mounts=[
                        lifesciences.Mount(disk="data", path="/data")
                    ]
                ),
                lifesciences.Action(
                    image_uri="google/cloud-sdk:slim",
                    commands=[
                        "gsutil", "-m", "cp", "-r", "/data/output/*", output_gcs_path
                    ],
                    mounts=[
                        lifesciences.Mount(disk="data", path="/data")
                    ]
                )
            ],
            resources=lifesciences.Resources(
                virtual_machine=lifesciences.VirtualMachine(
                    machine_type=machine_type,
                    boot_disk_size_gb=disk_size_gb,
                    service_account=lifesciences.ServiceAccount(
                        email=self.service_account or "default"
                    )
                )
            )
        )
        
        request = lifesciences.RunPipelineRequest(
            parent=self.parent,
            pipeline=pipeline
        )
        
        operation = self.lifesciences_client.run_pipeline(request=request)
        return operation.name
    
    def get_operation_status(self, operation_name: str) -> Dict[str, Any]:
        """Get status of a pipeline operation.
        
        Args:
            operation_name: Operation name.
        
        Returns:
            Operation status.
        """
        from google.longrunning import operations_pb2
        
        request = operations_pb2.GetOperationRequest(name=operation_name)
        operation = self.lifesciences_client.operations_client.get_operation(request)
        
        return {
            "name": operation.name,
            "done": operation.done,
            "error": str(operation.error) if operation.error.code else None,
            "metadata": dict(operation.metadata) if operation.metadata else {}
        }
    
    def wait_for_completion(
        self,
        operation_name: str,
        poll_interval: int = 60,
        timeout: int = 7200
    ) -> bool:
        """Wait for pipeline completion.
        
        Args:
            operation_name: Operation name.
            poll_interval: Seconds between checks.
            timeout: Maximum wait time.
        
        Returns:
            True if successful, False otherwise.
        """
        start_time = time.time()
        
        while True:
            status = self.get_operation_status(operation_name)
            
            if status["done"]:
                return status["error"] is None
            
            if time.time() - start_time > timeout:
                return False
            
            time.sleep(poll_interval)
    
    def download_results(
        self,
        gcs_output_path: str,
        local_output_dir: Path
    ) -> None:
        """Download results from GCS.
        
        Args:
            gcs_output_path: GCS path to results.
            local_output_dir: Local directory.
        """
        bucket_name, prefix = self._parse_gcs_path(gcs_output_path)
        bucket = self.storage_client.bucket(bucket_name)
        
        local_output_dir.mkdir(parents=True, exist_ok=True)
        
        blobs = bucket.list_blobs(prefix=prefix)
        
        for blob in blobs:
            relative_path = blob.name[len(prefix):].lstrip("/")
            local_path = local_output_dir / relative_path
            local_path.parent.mkdir(parents=True, exist_ok=True)
            
            blob.download_to_filename(str(local_path))
    
    def _parse_gcs_path(self, gcs_path: str) -> tuple:
        """Parse GCS path into bucket and prefix."""
        path = gcs_path.replace("gs://", "")
        parts = path.split("/", 1)
        bucket = parts[0]
        prefix = parts[1] if len(parts) > 1 else ""
        return bucket, prefix


def create_nextflow_config() -> str:
    """Generate Nextflow config for GCP."""
    return """
// Nextflow configuration for Google Cloud Life Sciences

profiles {
    gcp {
        process.executor = 'google-lifesciences'
        google.location = 'us-central1'
        google.region = 'us-central1'
        google.project = 'your-project-id'
        
        process {
            machineType = 'n1-standard-8'
            disk = '100 GB'
        }
        
        workDir = 'gs://your-bucket/work'
    }
}

process {
    withName: 'CNV_DETECTION' {
        cpus = 8
        memory = '32 GB'
        time = '6h'
    }
    
    withName: 'ANNOTATION' {
        cpus = 4
        memory = '16 GB'
        time = '2h'
    }
    
    withName: 'ML_PREDICTION' {
        cpus = 4
        memory = '16 GB'
        time = '1h'
    }
}

docker.enabled = true
"""


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run NRXN1 analysis on GCP")
    parser.add_argument("--project", required=True, help="GCP project ID")
    parser.add_argument("--input", required=True, help="GCS path to input VCF")
    parser.add_argument("--output", required=True, help="GCS path for output")
    parser.add_argument("--wait", action="store_true", help="Wait for completion")
    
    args = parser.parse_args()
    
    runner = GCPPipelineRunner(project_id=args.project)
    operation = runner.run_pipeline(args.input, args.output)
    print(f"Started pipeline: {operation}")
    
    if args.wait:
        print("Waiting for completion...")
        success = runner.wait_for_completion(operation)
        print(f"Pipeline {'succeeded' if success else 'failed'}")
