# Pistachio Growth Model Shiny App

An interactive Shiny application for visualizing **pistachio nut development using Growing Degree Days (GDD)**.  
This app implements a suite of statistically optimized growth models covering **kernel dry weight, kernel area, nut area, texture (hull/shell/kernel firmness), and nut + kernel color trajectories (L\*, a\*, b\*)**, derived from multi-year phenotyping data.

The application provides an intuitive interface for researchers, growers, and students to explore phenological trends, compare cultivars, and understand how thermal time influences nut development.

> **Modeling methods and formulas are summarized from the Pistachio Development Models document** .

---

## 🌱 Features

### **Core Functionality**
- Visualize predicted pistachio growth curves for multiple traits across GDD.
- Overlay observed field data for comparison.
- Switch between cultivars, sites, and trait categories.
- Explore **live** and **historical** GDD-based growth trajectories.
- Interactive hover labels for precise GDD and predicted trait values.

### **Color & Texture Visualization**
- Converts L\*a\*b\* predictions to **true RGB colors** with dynamic color boxes.
- Shows real-time color readout when hovering over plots.
- Displays kernel and nut firmness (texture) changes across the season.

### **Historical GDD Integration**
- Integrates 1980–2020 historical climate data (e.g., Cal-Adapt Livneh) to compute:
  - historical daily GDD  
  - historical cumulative GDD  
- Compares *current-year* nut development to historical averages.

---

## 📊 Model Summary (Methods)

The growth models were developed through a complete statistical workflow including **exploratory data analysis**, **cross-validation**, **transformations**, **mixed-effects modeling**, and **visual diagnostics** .

### **1. Model Selection Workflow**
- Analysis performed in **R 4.3.1**.
- Histograms and scatterplots used to assess raw data trends.
- **5-fold cross-validation** used to compare:
  - linear regression  
  - polynomial regressions (degrees 1–7)
- Best models chosen using:
  - **Lowest MSE**  
  - **Highest R²**  
  - **AIC/BIC** for transformed models  
- Residual and Q-Q plots guided the choice of:
  - Box-Cox transformations  
  - log transforms  
  - root transforms

### **2. Mixed-Effects Modeling**
A random effect for **year** was added to all polynomial models to account for repeated measurements and seasonal variability.

General model form:

\[
\hat{Y} = \beta_0 + \beta_1 x + \beta_2 x^2 + \cdots + \beta_j x^j
\]

where \(x\) = GDD (°F).

### **3. Back-Transformation**
After prediction, all transformed variables were **back-transformed into original biological units** for visualization.

---

## 🌰 Trait-Specific Models

### **Kernel Dry Weight**
- Fitted using a **cubic smoothing spline** due to polynomial overfitting.
- Hyperparameters:  
  - λ = 0.48021  
  - Degrees of Freedom = 3  
- Produces smooth, biologically realistic curves for kernel mass accumulation.

### **Kernel Area**
Second-degree polynomial with log-transformed predictor:

\[
\hat{Y}_{kernel\ area} =
-11001.714 + 2587.0618\log(x) - 149.9839\log(x)^2
\]

### **Nut Area**
Third-degree polynomial with log-transformed predictor:

\[
\hat{Y}_{nut\ area} =
-10917.8286 + 4438.07059\log(x)
-587.23897\log(x)^2 + 25.97156\log(x)^3
\]

### **Texture: Hull, Shell, Kernel**
- Kernel texture: Box-Cox transformed 2nd-degree polynomial (λ = −0.1414)  
- Shell texture: fourth-degree polynomial with log transformation  
- Hull texture: third-degree polynomial with log transformation  

These models capture firmness increases during Stages II → IV, as shown in the figures.

### **Nut Color (L\*, a\*, b\*)**
- **L\***: fifth-degree polynomial with log-transformation  
- **a\***: fourth-degree polynomial with log-transformation and shift  
- **b\***: fifth-degree Box-Cox model (λ = −0.26263)

### **Kernel Color**
- L\*: Box-Cox transformed 4-term polynomial (λ = −0.3030303)  
- a\*: third-degree log-transformed polynomial  
- b\*: third-degree log-transformed polynomial  

Color trajectories align with physiological color changes throughout fruit maturity.

---

## 🧭 App Layout

### **Live Predictions Tab**
- Select cultivar, trait, and measurement type.
- Display model predictions + observed points.
- Dynamic hover labels for precise values.

### **Historical Comparison Tab**
- Select location (lat/long) and date range.
- Compute historical vs current GDD curves.
- Overlay historical predicted growth trajectories.

### **Color Viewer**
- Converts predicted L\*, a\*, b\* values into RGB.
- Displays a color box that updates on hover to show real-time predicted nut or kernel color.

### **Model Details**
- Outlines formulas, assumptions, and transformations.

---
```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
