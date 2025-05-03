import shutil
import sys

def _check_dependencies(logger):
        """Check if all required external tools are available"""
        required_tools = ['gffread', 'bedtools', 'dustmasker', 'makeblastdb', 'tblastn', 'tfasty36']
        missing_tools = []
        
        for tool in required_tools:
            if shutil.which(tool) is None:
                missing_tools.append(tool)
        
        if missing_tools:
            logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            sys.exit(1)
        else:
            logger.info("All required tools are available")
