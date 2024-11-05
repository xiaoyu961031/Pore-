Welcome to the Pore+ Encoder repository, an enhanced feature engineering tool for metal-organic frameworks (MOFs) that builds upon traditional pore descriptors. The Pore+ encoder provides a refined set of descriptors designed to capture both the geometric and chemical heterogeneity of MOF pore spaces, enabling accurate and interpretable machine learning predictions for applications such as gas capture and separation.

Pore+ offers three key groups of descriptors: curated traditional geometric descriptors (Geo*), multi-probe-visited geometry (mpv-Geo), and advanced surface area metrics (SA). These descriptors aim to address the limitations of conventional pore features by incorporating more nuanced characterizations of the pore environment, including elementally and atom-type-specific contributions to surface area.

The Pore+ descriptors have been used successfully to predict selective adsorption of methyl iodide in MOFs, demonstrating improved model accuracy and interpretability compared to traditional descriptors. By providing a detailed representation of pore features, Pore+ helps researchers understand the critical factors driving adsorption performance and supports the rational design of new MOFs with optimized properties.

Elemental-decomposed srface area (e-SA):
To explore elemental decomposition, proceed with the example CD41/ provided by Prof. Linjiang Chen. For custom structures, supply a .mol file along with the appropriate probe and crystal parameters.

Atomic-type-decomposed surface area (at-SA):
For atomic-type decomposition, navigate to atomic-type/. If using your own database, adjust the element charge dictionary to align with user-defined requirements.

Multi-probe-visited porosity (mpv-Geo):
This feature is not included directly in the repository, as it can be easily derived from Zeo++ calculations. Ensure the results calculated via various probe sizes are processed prior to applying PCA for dimensionality reduction.

Special acknowledgement to Professor Linjiang Chen, from Precision and intelligent chemistry @ USTC & BHAM, for modifing code that decomposes surface area for POCs.
Special acknowledgement to Professor Tina Duren, from Univeresity of Bath, for providing the code prototype.

Referencee upcoming: The work is now under review for submission to Separation and Purification Technology. Xiaoyu Wu, Xianyu Song, Linjiang Chen*, Chunyi Yu, Liangdan Zhao, Mingrui Zuo, Chenrui Li, Heechae Choi, Jianwen Jiang* and Lifeng Ding*. Accurate and Interpretable Machine Learning with Pore+ Descriptors for Iodide Capture in Metal-Organic Frameworks
