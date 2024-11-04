FROM bioconductor/bioconductor_docker:3.19-R-4.4.0

# Add a non-root user and create the R library directory
RUN useradd -m cascade-user && \
    mkdir -p /home/cascade-user/Library/Frameworks/R.framework/Resources/site-library && \
    chown -R cascade-user:cascadeuser /home/cascade-user

# Set the R library path to the new directory
ENV R_LIBS_USER=/home/cascade-user/Library/Frameworks/R.framework/Resources/site-library

# Switch to the non-root user
USER cascade-user
WORKDIR /home/cascade-user

# Install R dependencies
RUN Rscript -e "devtools::install_github('adafede/cascade')"
# RUN Rscript -e "install.packages('cascade', repos = c('https://adafede.r-universe.dev', 'https://bioc.r-universe.dev', 'https://cran.r-universe.dev'))"

# Disable healthcheck (if you really want to disable it)
HEALTHCHECK NONE

# Define default command (commented out)
# CMD ["Rscript -e 'cascade::todo()'"]
