# pull ubuntu image
FROM julia:latest

# Install Julia package Wflow
RUN julia -e 'using Pkg; Pkg.add("Plots");' 

COPY . .

ENTRYPOINT ["julia"]
CMD ["./src/fontcreator.jl", "-p"]