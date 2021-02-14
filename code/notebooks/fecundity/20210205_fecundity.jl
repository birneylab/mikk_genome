### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 853a91fa-67b4-11eb-3647-9fe6553db874
begin
	using Pkg
	Pkg.add("XLSX")
end	

# ╔═╡ fce28616-67bb-11eb-1bf8-f71f5a8ffd4e
using Queryverse, XLSX, DataFrames, Plots

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
cd("/Users/brettell/Documents/Repositories/mikk_genome")

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
> * The day of egg collection during the week matters: the fish produce least on Tuesdays and Wednesdays; best on Thursday and Fridays; Monday are medium.
> * When females had no eggs, there's a \"0\". No entry means they weren't collected.
"

# ╔═╡ b6894624-67bf-11eb-2558-9b1b85cfd06d
readdir("data/fecundity")

# ╔═╡ 20097144-69e3-11eb-14f5-2bd787769258
df_quant = XLSX.readtable("data/fecundity/20210205_quantitative.xlsx", "Blatt1", header = true)

# ╔═╡ 359e1004-69e6-11eb-38b8-9174cf4b876d


# ╔═╡ 2d08df26-69e6-11eb-2a76-2716c7502d89


# ╔═╡ bb196aca-69e5-11eb-210e-6fe2a2290cf4


# ╔═╡ b1cebdce-69e5-11eb-3530-b9e754bd37cc


# ╔═╡ a3efa944-69e4-11eb-19c0-3d1687f0b1fa
df = DataFrame(df_quant, names = df_quant[1 :])

# ╔═╡ 58f06ea2-69e5-11eb-2237-3dc20d70f761


# ╔═╡ 07bce77c-69e5-11eb-064d-e7a1ff1ecf74
df_quant[1, :]

# ╔═╡ 3ebe1282-69e5-11eb-1f40-952f4d7ec10d


# ╔═╡ 3a1464fc-69e5-11eb-2d63-bf9417c64254


# ╔═╡ 0a82a8ac-69e5-11eb-3497-8d2052d4729b


# ╔═╡ f6f5466e-69e4-11eb-04c6-d794277f48d4


# ╔═╡ f47e7fcc-69e4-11eb-1dcd-8169a6ed8925


# ╔═╡ 9a0224b8-69e4-11eb-3f8d-9b61784d9ff0


# ╔═╡ 706ad370-69e4-11eb-3ea5-3d0ada9b7208


# ╔═╡ 12026668-69e4-11eb-06d2-a5b896f0d057


# ╔═╡ 4f39babe-67cc-11eb-3057-b51e3539912a


# ╔═╡ 1d60ec6a-67d9-11eb-29f8-4f2f3c2d34de


# ╔═╡ Cell order:
# ╟─ee9ea220-67b0-11eb-026f-752712eb67e0
# ╟─fbeea218-67b0-11eb-026d-e3857af31f42
# ╠═00a95b4a-67b1-11eb-0ee6-eb06d4d310d3
# ╠═7feec574-67b0-11eb-3d29-6d6dff03ea05
# ╠═7e92c610-67b4-11eb-10bd-c31b78949424
# ╠═853a91fa-67b4-11eb-3647-9fe6553db874
# ╠═fce28616-67bb-11eb-1bf8-f71f5a8ffd4e
# ╠═e12b13e8-67b1-11eb-04b2-5f5de4eab8a8
# ╟─f593f3e4-67c1-11eb-30ad-f52d6924078c
# ╠═b6894624-67bf-11eb-2558-9b1b85cfd06d
# ╠═20097144-69e3-11eb-14f5-2bd787769258
# ╠═359e1004-69e6-11eb-38b8-9174cf4b876d
# ╠═2d08df26-69e6-11eb-2a76-2716c7502d89
# ╠═bb196aca-69e5-11eb-210e-6fe2a2290cf4
# ╠═b1cebdce-69e5-11eb-3530-b9e754bd37cc
# ╠═a3efa944-69e4-11eb-19c0-3d1687f0b1fa
# ╠═58f06ea2-69e5-11eb-2237-3dc20d70f761
# ╠═07bce77c-69e5-11eb-064d-e7a1ff1ecf74
# ╠═3ebe1282-69e5-11eb-1f40-952f4d7ec10d
# ╠═3a1464fc-69e5-11eb-2d63-bf9417c64254
# ╠═0a82a8ac-69e5-11eb-3497-8d2052d4729b
# ╠═f6f5466e-69e4-11eb-04c6-d794277f48d4
# ╠═f47e7fcc-69e4-11eb-1dcd-8169a6ed8925
# ╠═9a0224b8-69e4-11eb-3f8d-9b61784d9ff0
# ╠═706ad370-69e4-11eb-3ea5-3d0ada9b7208
# ╠═12026668-69e4-11eb-06d2-a5b896f0d057
# ╟─4f39babe-67cc-11eb-3057-b51e3539912a
# ╟─1d60ec6a-67d9-11eb-29f8-4f2f3c2d34de
