= Preface

This book is written for people who use, or are planning to use, computational fluid dynamics (CFD) in their work, research, study or recreation time.

Its purpose is to provide the _knowledge_ to help the reader build their confidence to undertake CFD analysis, repeatedly, to a defined standard, in a timely manner#footnote[from Robin Hoyle's definition of 'competence', in _Complete training: from recruitment to retirement_, 2013, p69.]. It forms part of our mission to _make CFD accessible and inclusive_, alongside our work to manage, maintain and develop the free, open source software OpenFOAM on behalf of #link("https://openfoam.org")[The OpenFOAM Foundation Ltd].

The book taps into our experience in CFD, in particular with OpenFOAM. At the time of publication in 2022, it is 33 years since Henry Weller (HW) wrote the first lines of code of the “Field Operation and Manipulation” (FOAM) software which became OpenFOAM in 2004.

Since 2008 Chris Greenshields (CG) has provided CFD training with OpenFOAM, delivering over 300 courses (over 650 days) to 3,000 participants around the world. The training continues today at #link("https://cfd.direct/training")[CFD Direct Ltd], teaching repeatable _procedures_ and _workflows_ and the use of OpenFOAM and supporting software to deliver reliable CFD solutions.

This book does not replace training, but provides the supporting knowledge of the physics and numerical methods that underpin CFD. The user of any CFD software should benefit from reading this book, but they can only relate its contents back to the software if the code is open source, as it is with OpenFOAM.

CFD is challenging because it combines five sub-disciplines which are complex in their own right: (1) creating a computational mesh for the problem geometry; (2) fluid dynamics and physical modelling; (3) numerical methods; (4) data processing; and, (5) computers and programming.

This book concentrates mainly on fluid dynamics, modelling and numerical methods, finishing with a some examples that briefly describe some mesh generation and data processing. It presents:

- general equations of fluid dynamics that are the basis of a CFD solution;
- common physical models, including turbulence modelling;
- the finite volume method to generate equations in linearised matrix form;
- methods to compute solutions of the matrix equations;
- algorithms used to couple systems of equations;
- the practical application and numerical implementation of common boundary conditions;
- some sample problems.

We named the book _Notes on Computational Fluid Dynamics_ because we present the subject matter as a set of short topics, or _notes_, spanning exactly two pages of the print book, and one web page of this online version. Each topic appears in its entirety on the left and right pages of the opened, printed book. We use this format to break down the material into small pieces that are easier for the reader to digest.

The book presents CFD from the perspective of an engineer (or scientist), rather than a mathematician. CFD is an application of science so requires the qualities of an engineer to design, analyse, build and test, in order to create and run CFD simulations and interpret the results.

The finite volume method is the principal numerical method used in CFD. It is conceptually an engineering method that uses a computational mesh to define control volumes for physical processes. It is presented as such here, rather than as a mathematical construct which views the mesh as regions which specify distributions of variables.

Of course, we must use mathematics to describe physics, models and methods, but we do so in a way which is more familiar to engineers. Mathematics is just another form of language which we combine with verbal language in a manner that we hope conveys sufficient meaning while avoiding both ambiguity and unnecessary complexity.

For example, we describe a vector as a physical property, such as force and velocity, which has magnitude and direction, illustrated using an arrow. We do not describe it as a "three-dimensional vector space over the field of real numbers", using mathematical notation $in Re^3$.

We avoid the excessive use of sub- and superscripts within our mathematical notation. For example, when we describe an increase in time $t$ by an increment $Delta t$, we write $t := t + Delta t$, where "$:=$" means that the right hand side is evaluated and assigned to the variable on the left. This avoids the need to distinguish the "old" and "new" $t$ using an index $i$ by $t^(i+1) = t^i + Delta t$.

Numerical methods are ultimately programmed into a computer, so it is better that the mathematical notation represents typical programming language. The assignment notation "$t :=$", for example, clearly indicates a change to $t$ already stored in computer memory.

We have tried to reduce the need for the reader to search back through the book for information, to make it easier to follow. One way we hope to achieve this is by reminding the reader what the mathematical symbols represent, _e.g._ specific internal energy $e$ is defined in Sec. 2.15, but the reader is reminded what $e$ represents in Sec. 2.17 and Sec. 6.13.

When a topic refers back to something previously discussed, we have included references to direct the reader to the relevant section, page or equation number. To help the reader find something quickly, we have provided a fairly extensive index and list of symbols.

The book is designed to be self-contained so that the reader does not need to search through other texts to find some missing information. We have therefore dispensed with a long reference list at the end of the book in favour of footnotes within the pages themselves. Footnotes clarify minor points and present the origins of the science that is being used in CFD.

We use diagrams, graphs and images liberally to illustrate the topics we present. Each 2-page topic typically includes 1-2 figures, alongside equations. Some thematic figures are reused to reinforce certain concepts, redrawn with labelling relevant to the particular topic.

Each chapter includes one or more summary topics to give an overview of the main concepts presented within it.

Music has always been playing in the background during our work in CFD. As a light-hearted distraction, we chose to begin each chapter with relevant musical lyrics, which we think are just a poetic as anything found in literature.

The book is shaped by the people that taught, encouraged and supported us, past and present. CG would like to thank: Pat Leevers, Colin Wood and Bill Knox, for teaching him the value of writing, using diagrams and consistency; Richard Jones for many stimulating discussions about effective teaching; Steve Camp for his infectious enthusiasm and insightful questions; and, Akshai Runchal, Nils Basse and Bill Jones for their correspondence on some important topics.

HW would like to thank Chris Marooney, Bill Jones and David Gosman for their teaching and guidance in CFD.

Thanks to our colleagues Will Bainbridge and Jenya Collings for their support and, in Jenya's case, her assistance in publishing the book.

We hope you find the book useful.

Chris Greenshields and Henry Weller \
_CFD Direct Ltd_ / _The OpenFOAM Foundation Ltd_ \
April 2022.
