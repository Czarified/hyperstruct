This page is a blog post of sorts. It won't be around forever. It's just a way for Benjamin Crews to document
and potentially share his thoughts about all this.

# Why?

An immediate question, if you're not the sort to get excited about python code that automates structural analysis,
is "Why are you doing this...?" And the answer is pretty simple: I'm a huge aircraft nerd, and this is a passion
project of mine. I love advanced analytical concepts, dusting off old tombs of papers and documentation, and
making tools that can help people with advanced engineering topics. Additionally, I truly believe that the
aerospace industry would be better off if there were more open-source libraries and utilities tailored to our
unique requirements, written by experts, and hopefully validated against existing data or scientific research.
I don't think I'm an expert, and I'm definitely no scientist, but I'll do the best I can!

So, why did I choose _this_ exact implementation and tool to work on... is this really needed in the industry?
Well... things start to get more complicated here. The Society of Allied Weights Engineers (SAWE) publishes
many great papers, standards, procedures, and handbooks. I'm still very new to Weights Engineering (I became an
SAWE member in 2024), so I have a lot to learn here, but I hope through this project I can grow into an expert!
In the world of Weights Engineering, there are many procedures to estimate the weight of aircraft and aircraft
components. Most of these procedures, however, are rooted in historical data and many assumptions of design
practices. They often fail to answer questions about _how_ design changes affect the weight predictions.

## Basic Weight Approximations

As and example, consider the basic equations from Raymer. Raymer offers an Approximate Group Weights method,
which gives very simple formulas based on linear regression of historical aircraft data. All we need is the
wetted areas of geometry, and an idea of the aircraft's mission role. For a fighter aircraft this could work
out as follows:

$$
W_{wing} = 9.0 * S_{exposed planform} = 9.0 * 280.0 = 2520.0 [lbs]
W_{horz} = 4.0 * S_{exposed planform} = 4.0 * (2*32.2) = 257.6 [lbs]
W_{vert} = 4.3 * S_{exposed planform} = 4.3 * 31.85 = 137.0 [lbs]
W_{fuse} = 4.8 * S_{wetted area} = 4.8 * 128.7 = 617.8 [lbs]
W_{gear} = 0.033 * TOGW = 0.033 * 37500 = 1237.5 [lbs]
W_{inst.eng} = 1.3 * W_eng = 1.3 * 4100.0 = 5330.0 [lbs]
W_{misc} = 0.17 * TOGW = 0.17 * 37500 = 6375.0 [lbs]

OEW = 16,475 [lbs]
$$

So that's a bit of math, but in the end we end up with an answer _kinda_ close to the empty weight of an F-16
Block 50 (18,900 lbs). Clearly it's off, but it's also incredibly simple and straight-forward. This is great
for the early conceptual design phase, and especially when the structural engineers aren't yet involved.

Stepping up in fidelity Raymer proposes a statistical approach with many more factors and inputs. This
method is theoretically closer to the true value, but again assumes typical construction methods and standard
design practices. In the evolving world of weight-critical structures and searches for efficiency, this becomes
much less useful, much earlier in the design process. As an example, take the F-16 Horizontal Stabilizer: A
design study was conducted in 1982 on the weight impacts of (at the time) advanced composite materials
{cite}`buther_non-honeycomb_1982`. This study shows that a non-honeycomb metallic stabilizer weighed
approximately 167 lbs (compare that to half our previous estimate for a single stabilizer), but with advanced
materials it could be reduced to 158 lbs!

To my knowledge, most weight approximation methods follow a similar procedure, trading off increasing
mathematical and statistical complexity for increasingly accurate estimates. Statistical methods have factors
for things like advanced materials, or configuration changes. This is one of the duty of the Weights Engineer
during the Design Process.

## Digging Deeper

At a certain point in the design cycle, structural designers start asking questions like "Is a stringer or longeron
system more efficient? How close do these frames _really_ need to be? What's the weight tradeoff of
this wing being mid-fuselage mounted vs bottom-mounted? How much weight gets added to the fuselage when we add an
active sweep mechanism into the wing?" These questions are all quite complicated to answer
and often require a (team of) Stress/Structural Engineer(s) to conduct detailed analysis studies over a much longer
period of time than a simple factor lookup and summation of line items. And yet, these can be very important
questions to the overall design! Often, teams can rely upon the expertise of those that came before, but
that's not always correct and not always available.

And this scene is where SWEEP comes in.

> ... the risk is high of adoptin a less than optimum basic configuration at the very inception of
> a program,adversely affecting system effectiveness throughout the life of the design.

SWEEP's documentation is extensive. I haven't even read through it all. The executive summary, however, is
the most useful at, well... _summarizing_ what it does and why {cite}`vol1_2850`. In general, structural
weight is a **huge** driving function of total vehicle weight, and therefore effectivity. To use existing
methods, you _need_ derived statistical parameter, relevant to your mission and flight profiles. What if
you start pushing the limits beyond the data collected? Fly higher, fly faster, pull harder g's, or even
fly into the realm where thermomechanical effects begin to drive the structural sizing. You need an analytical
solution instead of statistical, or perhaps you could iteratively build planes on the fringe of your dataset
and proceed to collect more data each time.

Additionally, in modern times AI/ML and Large Language Models are
commonly used in place of statistical derivations and factors. These have the same problem. To build something
beyond the existing dataset, you inherently assume that the trends continue. Maybe they do, but a physics-based
approach will give more confidence in the estimates.

## ... but SWEEP already exists...?

Maybe! I certainly don't have access to it. Given it was written in 1974, in FORTRAN, for some specific
mainframe hardware systems, it's unlikely that even with the original source code we could execute it.
Updating this into a modern language opens up the tools and methods for today's, and the future's, aircraft
designers and engineers! So far, I haven't been able to find a legible copy of the original source code.
With the original code, perhaps AI could convert it into a full executable library, but I have my doubts.
So, this leaves us with a proverbial mountain of method documentation, and very long list of work to cobble
pieces off that mountain and build into some functional structure...

I selected python as the language for this project, because I'm fluent and love it. If python's execution
speed becomes an issue, such that languages like C# or Julia need to be considered, this project is already
a huge success. Until that time, python should be fast enough and have enough features to executed what I need.

SWEEP is documented and archived by the Defence Technical Information Center.
Search for the relevant technical number, or maybe the links below will help...
{cite}`vol1_2850`
{cite}`vol2_2851`
{cite}`vol2_2852`
{cite}`vol2_2853`
{cite}`vol3_2854`
{cite}`vol3_2855`
{cite}`vol4_2856`
{cite}`vol5_2857`
{cite}`vol5_2858`
{cite}`vol6_2859`
{cite}`vol6_2860`
{cite}`vol6_2861`
{cite}`vol6_2862`
{cite}`vol6_2863`
{cite}`vol6_2864`
{cite}`vol6_2865`
{cite}`vol6_2866`
{cite}`vol7_2867`
{cite}`vol7_2868`
{cite}`vol8_2869`
{cite}`vol9_2871`
{cite}`vol11_2873`

# References

```{bibliography}

```
