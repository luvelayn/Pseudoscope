import logging
import os

def _setup_logging(output_dir):
        """Set up logging configuration"""
        logger = logging.getLogger('pseudoscope')
        logger.setLevel(logging.INFO)
        
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S')
        
        # Console handler
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        
        # File handler
        fh = logging.FileHandler(os.path.join(output_dir, 'pseudoscope.log'))
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        
        return logger
