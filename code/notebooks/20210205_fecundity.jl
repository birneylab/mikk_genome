### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 853a91fa-67b4-11eb-3647-9fe6553db874
begin
	using Pkg
	Pkg.add("DataVoyager")
end	

# ╔═╡ fce28616-67bb-11eb-1bf8-f71f5a8ffd4e
using Queryverse, DataFrames, Plots

# ╔═╡ 4f39babe-67cc-11eb-3057-b51e3539912a
using DataVoyager

# ╔═╡ ee9ea220-67b0-11eb-026f-752712eb67e0
md"
# Fecundity plots
"

# ╔═╡ fbeea218-67b0-11eb-026d-e3857af31f42
md"
## Setup
"

# ╔═╡ 00a95b4a-67b1-11eb-0ee6-eb06d4d310d3
md"
### Set working directory
"

# ╔═╡ 7feec574-67b0-11eb-3d29-6d6dff03ea05
cd("/hps/research1/birney/users/ian/mikk_paper/mikk_genome")

# ╔═╡ 7e92c610-67b4-11eb-10bd-c31b78949424
md"
### Load packages
"

# ╔═╡ e12b13e8-67b1-11eb-04b2-5f5de4eab8a8
md"
### Read in data
"

# ╔═╡ f593f3e4-67c1-11eb-30ad-f52d6924078c
md"
#### Info on data

____________________________

**Semiquantitative**:
* Document `MIKK active lines 9-20 fecundity.xlsx` sent from FL to TF on 20200911.
* Saved as `data/fecundity/20210205_semiquantitative.xlsx`
> * \"Average\" over several weeks of egg collection and in most cases several females per strain.
> * Twice, once in Feb/Mar 2019 then again summer 2020. 
> * In both cases the fishes were about 6 months old, at the peak of their egg production.
> * Not a well-controlled experiment; not enough space to set up tanks with fixed numbers of fish to go into real quantification. 
> * Therefore only a rough estimate of fecundity.

____________________________

**Quantitative**:

* Document `MIKK active lines F15-16 females.xlsx` sent from FL to TF on 20200918.
* Saved as `data/fecundity/20210205_quantitative.xlsx`
> * Shows numbers of eggs collected on a given day from a number of females. 
> * Counted eggs from strains at F15 and F16 (i.e. eggs are F16 and F17 respectively). 
> * Average of the fishes is indicated in the table. F15 fish were about 6 months older than F16 eggs. 
> * The day of egg collection during the week matters: the fish produce least on Tuesdays and Wednesdays; best on Thursday and Fridays; 
"

# ╔═╡ 23fc39e4-67c2-11eb-2b28-cbfd3515c5fb


# ╔═╡ b6894624-67bf-11eb-2558-9b1b85cfd06d
readdir("data/fecundity")

# ╔═╡ 2ac85402-67b2-11eb-0fed-7389ff6d7f85
XLSX.readxlsx("data/fecundity")

# ╔═╡ 1d60ec6a-67d9-11eb-29f8-4f2f3c2d34de
v = Voyager()

# ╔═╡ Cell order:
# ╟─ee9ea220-67b0-11eb-026f-752712eb67e0
# ╟─fbeea218-67b0-11eb-026d-e3857af31f42
# ╠═00a95b4a-67b1-11eb-0ee6-eb06d4d310d3
# ╠═7feec574-67b0-11eb-3d29-6d6dff03ea05
# ╠═7e92c610-67b4-11eb-10bd-c31b78949424
# ╠═853a91fa-67b4-11eb-3647-9fe6553db874
# ╠═fce28616-67bb-11eb-1bf8-f71f5a8ffd4e
# ╠═e12b13e8-67b1-11eb-04b2-5f5de4eab8a8
# ╠═f593f3e4-67c1-11eb-30ad-f52d6924078c
# ╠═23fc39e4-67c2-11eb-2b28-cbfd3515c5fb
# ╠═b6894624-67bf-11eb-2558-9b1b85cfd06d
# ╠═2ac85402-67b2-11eb-0fed-7389ff6d7f85
# ╠═4f39babe-67cc-11eb-3057-b51e3539912a
# ╠═1d60ec6a-67d9-11eb-29f8-4f2f3c2d34de
