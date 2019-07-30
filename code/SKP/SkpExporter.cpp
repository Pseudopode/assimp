#ifndef ASSIMP_BUILD_NO_EXPORT
#ifndef ASSIMP_BUILD_NO_SKP_EXPORTER

#include "SkpExporter.h"
#include <assimp/Bitmap.h>
#include <assimp/fast_atof.h>
#include <assimp/SceneCombiner.h>
#include <assimp/StringUtils.h>
#include <assimp/DefaultIOSystem.h>
#include <assimp/IOSystem.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>

#include <assimp/Exceptional.h>

#include <SketchUpAPI/sketchup.h>

#include <cassert>
#include <ctime>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

using namespace Assimp;

namespace Assimp {
namespace {

#define SU(api_function_call) {\
SUResult su_api_result = api_function_call;\
assert(SU_ERROR_NONE == su_api_result);\
}\

std::string GetString(const SUStringRef& string) {
    size_t length = 0;
    SU(SUStringGetUTF8Length(string, &length));
    std::vector<char> buffer(length + 1);
    size_t out_length = 0;
    SU(SUStringGetUTF8(string, length + 1, buffer.data(), &out_length));
    assert(out_length == length);
    return std::string(begin(buffer), end(buffer));
}

const double MeterToInch = 39.37007874015748;
const double MilliMeterToInch = 0.03937007874015748;


SUVector3D ComputeNormal(const aiMesh* mesh, const aiFace& face) {
    // Copied from GenFaceNormalsProcess.
    const aiVector3D* pV1 = &mesh->mVertices[face.mIndices[0]];
    const aiVector3D* pV2 = &mesh->mVertices[face.mIndices[1]];
    const aiVector3D* pV3 = &mesh->mVertices[face.mIndices[face.mNumIndices - 1]];
    const aiVector3D vNor = ((*pV2 - *pV1) ^ (*pV3 - *pV1)).NormalizeSafe();
    return { vNor.x, vNor.y, vNor.z };
}

SUMaterialRef LoadMaterial(aiMaterial* ai_material, SUModelRef model) {
    SUMaterialRef material = SU_INVALID;
    SU(SUMaterialCreate(&material));

    std::cout << "Material: " << ai_material->GetName().C_Str() << "\n";
    SU(SUMaterialSetName(material, ai_material->GetName().C_Str())); // TODO: Ensure is unique?

    aiColor4D ai_color;
    if (aiReturn_SUCCESS == ai_material->Get(AI_MATKEY_COLOR_DIFFUSE, ai_color)) {
        SUByte r = static_cast<SUByte>(std::round(255.0 * ai_color.r));
        SUByte g = static_cast<SUByte>(std::round(255.0 * ai_color.g));
        SUByte b = static_cast<SUByte>(std::round(255.0 * ai_color.b));
        SUByte a = static_cast<SUByte>(std::round(255.0 * ai_color.a));
        std::cout << "  aiColor: " << ai_color.r << ", " << ai_color.g << ", " << ai_color.b << ", " << ai_color.a << "\n";
        std::cout << "  SUColor: " << unsigned(r) << ", " << unsigned(g) << ", " << unsigned(b) << ", " << unsigned(a) << "\n";

        SUColor color{ r, g, b, a };
        SU(SUMaterialSetColor(material, &color));
    }

    ai_real opacity;
    if (aiReturn_SUCCESS == ai_material->Get(AI_MATKEY_OPACITY, opacity)) {
        // TODO: SUMaterialSetUseOpacity ?
        SU(SUMaterialSetOpacity(material, opacity));
    }

    if (ai_material->GetTextureCount(aiTextureType_DIFFUSE) > 0) {
        aiString ai_path;
        aiTextureMapping ai_mapping; // Handle only aiTextureMapping_UV
        aiTextureMapMode ai_mode; // Handle only aiTextureMapMode_Wrap
        auto result = ai_material->GetTexture(aiTextureType_DIFFUSE, 0,
            &ai_path, &ai_mapping, nullptr, nullptr, nullptr, &ai_mode);
        assert(result == aiReturn_SUCCESS);
        assert(ai_mapping == aiTextureMapping_UV);
        assert(ai_mode == aiTextureMapMode_Wrap);
        std::cout << "  Texture: " << ai_path.C_Str() << "\n";

        double scale_s = 1.0;
        double scale_t = 1.0;
        SUTextureRef texture = SU_INVALID;
        SU(SUTextureCreateFromFile(&texture, ai_path.C_Str(), scale_s, scale_t)); // TODO: Handle missing file.
        SU(SUMaterialSetTexture(material, texture));
    }

    //SU(SUModelAddMaterials(model, 1, &material));
    auto result = SUModelAddMaterials(model, 1, &material);
    assert(result == SU_ERROR_NONE);
    return material;
}

void SetModelOptions(SUModelRef model) {
    SUOptionsManagerRef options = SU_INVALID;
    SU(SUModelGetOptionsManager(model, &options));

    SUOptionsProviderRef unit_options = SU_INVALID;
    SU(SUOptionsManagerGetOptionsProviderByName(options,
        "UnitsOptions", &unit_options));

    SUTypedValueRef value = SU_INVALID;
    SU(SUTypedValueCreate(&value));

    // Length Units

    SU(SUTypedValueSetInt32(value, SU_LFORMAT_DECIMAL));
    SU(SUOptionsProviderSetValue(unit_options, "LengthFormat", value));

    SU(SUTypedValueSetInt32(value, SU_LUNIT_METER));
    SU(SUOptionsProviderSetValue(unit_options, "LengthUnit", value));

    SU(SUTypedValueSetInt32(value, 2)); // 0.00
    SU(SUOptionsProviderSetValue(unit_options, "LengthPrecision", value));

    SU(SUTypedValueSetBool(value, true));
    SU(SUOptionsProviderSetValue(unit_options, "LengthSnapEnabled", value));

    SU(SUTypedValueSetDouble(value, 10.0 * MilliMeterToInch));
    SU(SUOptionsProviderSetValue(unit_options, "LengthSnapLength", value));

    // Area Units
    
    SU(SUTypedValueSetInt32(value, SU_AUNIT_SQUARE_METER));
    SU(SUOptionsProviderSetValue(unit_options, "AreaUnit", value));

    // Volume Units

    SU(SUTypedValueSetInt32(value, SU_VUNIT_CUBIC_METER));
    SU(SUOptionsProviderSetValue(unit_options, "VolumeUnit", value));

    // Angle Units

    SU(SUTypedValueSetBool(value, true));
    SU(SUOptionsProviderSetValue(unit_options, "AngleSnapEnabled", value));

    SU(SUTypedValueSetInt32(value, 1));
    SU(SUOptionsProviderSetValue(unit_options, "AnglePrecision", value));

    SU(SUTypedValueSetDouble(value, 15.0));
    SU(SUOptionsProviderSetValue(unit_options, "SnapAngle", value));


    SURenderingOptionsRef rendering_options = SU_INVALID;
    SU(SUModelGetRenderingOptions(model, &rendering_options));

    // TODO: Turn off for production.
    SU(SUTypedValueSetBool(value, true));
    SURenderingOptionsSetValue(rendering_options, "DisplayInstanceAxes", value);

    SU(SUTypedValueRelease(&value));
}

void SetModelCameraIsoPerspective(SUModelRef model) {
    SUCameraRef camera = SU_INVALID;
    SU(SUCameraCreate(&camera));
    SUPoint3D eye{ -3.0 * MeterToInch, -3.0 * MeterToInch, 2.0 * MeterToInch };
    SUPoint3D target{ 0.0, 0.0, 0.5 * MeterToInch };
    SUVector3D up{ 0.0, 0.0, 1.0 };
    SU(SUCameraSetOrientation(camera, &eye, &target, &up));
    SU(SUCameraSetPerspectiveFrustumFOV(camera, 35)); // degrees
    SU(SUModelSetCamera(model, &camera));
}

void SetRenderingOptions(SURenderingOptionsRef rendering_options,
        const char* key, bool value) {
    SUTypedValueRef typed_value = SU_INVALID;
    SU(SUTypedValueCreate(&typed_value));
    SU(SUTypedValueSetBool(typed_value, value));
    SU(SURenderingOptionsSetValue(rendering_options, key, typed_value));
    SU(SUTypedValueRelease(&typed_value));
}

void SetRenderingOptions(SURenderingOptionsRef rendering_options,
        const char* key, int32_t value) {
    SUTypedValueRef typed_value = SU_INVALID;
    SU(SUTypedValueCreate(&typed_value));
    SU(SUTypedValueSetInt32(typed_value, value));
    SU(SURenderingOptionsSetValue(rendering_options, key, typed_value));
    SU(SUTypedValueRelease(&typed_value));
}

void SetRenderingOptions(SURenderingOptionsRef rendering_options,
        const char* key, SUColor value) {
    SUTypedValueRef typed_value = SU_INVALID;
    SU(SUTypedValueCreate(&typed_value));
    SU(SUTypedValueSetColor(typed_value, &value));
    SU(SURenderingOptionsSetValue(rendering_options, key, typed_value));
    SU(SUTypedValueRelease(&typed_value));
}

void SetModelStyle(SUModelRef model) {
    // TODO: The SketchUp C API doesn't support all style properties :(
    // Set the properties via RenderingOptions even though that will make the
    // active style marked as changed.

    // This is mostly replicating the "Simple Style" from SketchUp.

    SURenderingOptionsRef ro = SU_INVALID;
    SU(SUModelGetRenderingOptions(model, &ro));

    SetRenderingOptions(ro, "BackgroundColor", { 153, 189, 144, 255 });
    SetRenderingOptions(ro, "BandColor", { 0, 0, 0, 127 });
    SetRenderingOptions(ro, "ConstructionColor", { 0, 0, 0, 255 });
    SetRenderingOptions(ro, "FaceBackColor", { 164, 178, 187, 255 });
    SetRenderingOptions(ro, "FaceFrontColor", { 255, 255, 255, 255 });
    SetRenderingOptions(ro, "FogColor", { 153, 189, 144, 255 });
    SetRenderingOptions(ro, "ForegroundColor", { 0, 0, 0, 255 });
    SetRenderingOptions(ro, "GroundColor", { 210, 208, 185, 255 });
    SetRenderingOptions(ro, "HighlightColor", { 0, 1, 255, 255 });
    SetRenderingOptions(ro, "HorizonColor", { 0, 0, 0, 0 });
    SetRenderingOptions(ro, "LockedColor", { 255, 0, 0, 255 });
    SetRenderingOptions(ro, "SectionActiveColor", { 255, 135, 0, 255 });
    SetRenderingOptions(ro, "SectionDefaultCutColor", { 0, 0, 0, 255 });
    SetRenderingOptions(ro, "SectionDefaultFillColor", { 0, 0, 0, 255 });
    SetRenderingOptions(ro, "SectionInactiveColor", { 112, 105, 97, 255 });
    SetRenderingOptions(ro, "SkyColor", { 104, 195, 251, 255 });

    SetRenderingOptions(ro, "DrawHorizon", true);
    SetRenderingOptions(ro, "SilhouetteWidth", 2);
}

SUTransformation ColladaToSkpCoordinates() {
    SUTransformation scale_transformation;
    SU(SUTransformationScale(&scale_transformation, MeterToInch));

    SUPoint3D origin{ 0.0, 0.0, 0.0 };
    SUVector3D x_axis{ -1.0, 0.0, 0.0 }; // TODO: Is negative X-Axis correct?
    SUVector3D y_axis{ 0.0, 0.0, 1.0 };
    SUVector3D z_axis{ 0.0, 1.0, 0.0 };
    SUTransformation axes_transformation;
    SU(SUTransformationSetFromPointAndAxes(&axes_transformation,
        &origin, &x_axis, &y_axis, &z_axis));

    SUTransformation transformation;
    SU(SUTransformationMultiply(&axes_transformation, &scale_transformation,
        &transformation));

    return transformation;
}


class SkpExporter {
public:
    SkpExporter();

    void Export(const char* filepath, const aiScene* scene);

private:
    void LoadMaterials(const aiScene* scene);
    void LoadNodeMeshes(const aiNode* node, SUEntitiesRef entities);
    void LoadNodes(const aiScene* scene);

    void MeshToGeometryInput(aiMesh* mesh, SUGeometryInputRef input);
    void NodeToInstance(aiNode* node, SUEntitiesRef entities);

    // Collada to SKP coordinate system transformation.
    const SUTransformation transformation_{};

    const aiScene* scene_ = nullptr;
    SUModelRef model_ = SU_INVALID;
    std::vector<SUMaterialRef> materials_;
    std::vector<SUComponentDefinitionRef> definitions_;
};


SkpExporter::SkpExporter() : transformation_(ColladaToSkpCoordinates()) {
}

void SkpExporter::Export(const char* filepath, const aiScene* scene) {
    assert(scene);
    scene_ = scene;
    assert(SUIsInvalid(model_));
    model_ = SU_INVALID;
    SU(SUModelCreate(&model_));
    SetModelOptions(model_);
    SetModelStyle(model_);
    SetModelCameraIsoPerspective(model_);
    LoadMaterials(scene);
    LoadNodes(scene);
    SU(SUModelMergeCoplanarFaces(model_));
    SU(SUModelSaveToFile(model_, filepath)); // TODO: might be locked.
    SU(SUModelRelease(&model_));
}

void SkpExporter::LoadMaterials(const aiScene* scene) {
    for (size_t i = 0; i < scene->mNumMaterials; i++) {
        auto ai_material = scene->mMaterials[i];
        auto material = LoadMaterial(ai_material, model_);
        materials_.emplace_back(material);
    }
}

void SkpExporter::LoadNodes(const aiScene* scene) {
    SUEntitiesRef entities = SU_INVALID;
    SU(SUModelGetEntities(model_, &entities));
    NodeToInstance(scene->mRootNode, entities);
}

void SkpExporter::LoadNodeMeshes(const aiNode* node, SUEntitiesRef entities) {
    if (node->mNumMeshes == 0) return;

    SUGeometryInputRef input = SU_INVALID;
    SU(SUGeometryInputCreate(&input));

    for (size_t i = 0; i < node->mNumMeshes; i++) {
        auto mesh_index = node->mMeshes[i];
        aiMesh* mesh = scene_->mMeshes[mesh_index];
        MeshToGeometryInput(mesh, input);
    }
    SU(SUEntitiesFill(entities, input, true));

    SU(SUGeometryInputRelease(&input));

    // Fix the smooth bug in SULoopInputEdgeSetSmooth.
    size_t num_edges = 0;
    SU(SUEntitiesGetNumEdges(entities, false, &num_edges));
    std::vector<SUEdgeRef> edges(num_edges, SU_INVALID);
    SU(SUEntitiesGetEdges(entities, false, num_edges, edges.data(), &num_edges));
    for (auto& edge : edges) {
        bool is_soft = false;
        SU(SUEdgeGetSoft(edge, &is_soft));
        if (is_soft) {
            SU(SUEdgeSetSmooth(edge, true));
        }
    }
}

void SkpExporter::MeshToGeometryInput(aiMesh* mesh, SUGeometryInputRef input) {
    assert(mesh->HasFaces());
    assert(mesh->HasPositions());
    assert(mesh->HasNormals());

    std::cout << "Mesh: " << mesh->mName.C_Str() << "\n";

    size_t vertices_offset = 0, faces_offset = 0, edges_offset = 0, curves_offset = 0, arcs_offset = 0;
    SU(SUGeometryInputGetCounts(input, &vertices_offset, &faces_offset, &edges_offset, &curves_offset, &arcs_offset));

    std::cout << "  Vertex offset: " << vertices_offset << "\n";

    SUMaterialRef material = materials_.at(mesh->mMaterialIndex);

    for (size_t i = 0; i < mesh->mNumVertices; i++) {
        auto vertex = mesh->mVertices[i];
        SUPoint3D position{ vertex.x, vertex.y, vertex.z };
        SU(SUPoint3DTransform(&transformation_, &position));
        SU(SUGeometryInputAddVertex(input, &position));
    }

    for (size_t i = 0; i < mesh->mNumFaces; i++) {
        auto face = mesh->mFaces[i];
        //std::cout << "  Face: " << i << "\n";

        SULoopInputRef loop = SU_INVALID;
        SU(SULoopInputCreate(&loop));

        for (size_t j = 0; j < face.mNumIndices; j++) {
            size_t vertex_index = vertices_offset + face.mIndices[j];
            SU(SULoopInputAddVertexIndex(loop, vertex_index));
        }

        // Set edges to soft+smooth if vertex normals isn't in the same
        // direction as the face normal. This is the best SketchUp can do since
        // it doesn't allow explicit control over vertex normals.
        SUVector3D face_normal = ComputeNormal(mesh, face);
        //std::cout << "    Normal: " << face_normal.x << ", " << face_normal.y << ", " << face_normal.z << "\n";
        for (size_t j = 0; j < face.mNumIndices; j++) {
            // Get the index of the next vertex in the loop.
            size_t k = (j + 1) % face.mNumIndices;
            // Fetch the vertex normals.
            size_t ni1 = face.mIndices[j];
            size_t ni2 = face.mIndices[k];
            auto n1 = mesh->mNormals[ni1];
            auto n2 = mesh->mNormals[ni2];
            SUVector3D normal1{ n1.x, n1.y, n1.z };
            SUVector3D normal2{ n2.x, n2.y, n2.z };
            // Check if they are in the same direction as the face.
            bool same1 = false, same2 = false;
            SU(SUVector3DIsSameDirectionAs(&face_normal, &normal1, &same1));
            SU(SUVector3DIsSameDirectionAs(&face_normal, &normal2, &same2));

            //std::cout << "      V1 Normal: " << normal1.x << ", " << normal1.y << ", " << normal1.z << " same: " << same1 << "\n";
            //std::cout << "      V2 Normal: " << normal2.x << ", " << normal2.y << ", " << normal2.z << " same: " << same2 << "\n";

            // If either vertex normal is not in the same direction then assume
            // this edge should be soft+smooth.
            if (!same1 || !same2) {
                //std::cout << "    Soft+Smooth edge: " << j << "\n";
                size_t edge_index = j;
                SU(SULoopInputEdgeSetSoft(loop, edge_index, true));
                // TODO: SketchUp Bug! The smooth propery appear to be applied
                //   to the whole loop instead of just the given edge index.
                //SU(SULoopInputEdgeSetSmooth(loop, edge_index, true));
            }
        }

        size_t face_index = 0;
        SU(SUGeometryInputAddFace(input, &loop, &face_index));

        size_t num_uv_coords = 0;
        SUMaterialInput material_input{ num_uv_coords, {}, {}, material };
        SU(SUGeometryInputFaceSetFrontMaterial(input, face_index, &material_input));
    }
}

void SkpExporter::NodeToInstance(aiNode* node, SUEntitiesRef entities) {
    std::cout << "Node: " << node->mName.C_Str() << "\n";

    if (node->mNumChildren == 0 && node->mNumMeshes == 0) {
        return;
    }

    std::cout << "  Meshes: " << node->mNumMeshes << "\n";
    LoadNodeMeshes(node, entities);

    for (size_t i = 0; i < node->mNumChildren; i++) {
        auto child = node->mChildren[i];

        auto tr = child->mTransformation;
        //std::cout << "  Transformation:\n";
        //std::cout << tr.a1 << ", " << tr.a2 << ", " << tr.a3 << ", " << tr.a4 << "\n";
        //std::cout << tr.b1 << ", " << tr.b2 << ", " << tr.b3 << ", " << tr.b4 << "\n";
        //std::cout << tr.c1 << ", " << tr.c2 << ", " << tr.c3 << ", " << tr.c4 << "\n";
        //std::cout << tr.d1 << ", " << tr.d2 << ", " << tr.d3 << ", " << tr.d4 << "\n";

        // assimp is row-major while SketchUp is column-major
        const double scale = MeterToInch;
        SUTransformation transformation{
            tr.a1, tr.b1, tr.c1, tr.d1,
            tr.a2, tr.b2, tr.c2, tr.d2,
            tr.a3, tr.b3, tr.c3, tr.d3,
            // Beware that the Y and Z position is swapped to account for the YZ
            // axis difference.
            tr.a4 * scale, tr.c4 * scale, tr.b4 * scale, tr.d4,
        };

        SUGroupRef group = SU_INVALID;
        SU(SUGroupCreate(&group));
        SU(SUGroupSetName(group, child->mName.C_Str()));
        SU(SUGroupSetTransform(group, &transformation));
        SU(SUEntitiesAddGroup(entities, group));

        SUEntitiesRef child_entities = SU_INVALID;
        SU(SUGroupGetEntities(group, &child_entities));

        NodeToInstance(child, child_entities);
    }
}

} // namespace


// ------------------------------------------------------------------------------------------------
// Worker function for exporting a scene to Collada. Prototyped and registered in Exporter.cpp
void ExportSceneSKP(const char* pFile, IOSystem* pIOSystem,
        const aiScene* pScene, const ExportProperties* /*pProperties*/) {

    // <debug>
    std::cout << "\n";
    std::cout << "Metadata:\n";
    for (size_t i = 0; i < pScene->mMetaData->mNumProperties; i++) {
        auto key = pScene->mMetaData->mKeys[i];
        //auto value = pScene->mMetaData->mValues[i];
        std::cout << key.C_Str() << "\n";
    }
    std::cout << "\n";
    // </debug>

    SUInitialize();
    {
        SkpExporter exporter;
        exporter.Export(pFile, pScene);
    }
    SUTerminate();
}

} // end of namespace Assimp

#endif
#endif
